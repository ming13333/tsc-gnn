#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step L3-upgrade: 用公开人卒中 bulk (GSE16561, 外周血, 39 卒中 vs 24 对照)
跑 ssGSEA / AUCell 风格的模块活性打分，把"结构保守"升级为"表达激活保守"。

流程:
  1. 下载 GSE16561 series matrix (ILMN 探针 x 样本, RMA 归一化) + GPL6883 注释
  2. 探针 -> 基因符号, 折叠为基因 x 样本
  3. 从 human DoRothEA 提取 SOX10/CEBPB/GATA2 人源靶基因集
  4. 每样本每模块 ssGSEA 风格打分 (rank-mass fraction, [0,1])
  5. stroke vs control: Wilcoxon + Cliff's delta + BH 校正
"""
import urllib.request, gzip, os, json, re
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
BULK_CACHE = os.path.join(HERE, "gse16561_bulk.npz")
ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6883/annot/GPL6883.annot.gz"
MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16561/matrix/GSE16561_series_matrix.txt.gz"
DOROTHEA = os.path.join(os.path.dirname(HERE), "data", "dorothea", "human_dorothea_regulon.tsv")
TFS = ["SOX10", "CEBPB", "GATA2"]
NC_TFS = ["OLIG2", "PAX6", "NKX2-2"]  # 脑特异 TF，预期在血液中不应激活 -> 负对照


def download(url, timeout=180):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=timeout).read()


def parse_matrix(raw):
    txt = gzip.decompress(raw).decode("utf-8", "ignore")
    lines = txt.splitlines()
    table_begin = None
    title = desc = None
    for i, l in enumerate(lines):
        if l.startswith("!series_matrix_table_begin"):
            table_begin = i
        if l.startswith("!Sample_title\t"):
            title = l.split("\t")[1:]
        if l.startswith("!Sample_description\t"):
            desc = l.split("\t")[1:]
    hdr = lines[table_begin + 1].split("\t")
    samples = hdr[1:]
    probes, mat = [], []
    for l in lines[table_begin + 2:]:
        if l.startswith("!series_matrix_table_end"):
            break
        p = l.split("\t")
        probes.append(p[0].strip().strip('"'))
        vals = []
        for x in p[1:]:
            x = x.strip()
            if x in ("", "NA", "null"):
                vals.append(np.nan)
            else:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(np.nan)
        mat.append(vals)
    M = np.array(mat, dtype=float)  # probes x samples
    # grouping from title suffix
    groups = []
    for t in title:
        t = t.strip().strip('"')
        groups.append("Stroke" if t.endswith("Stroke") else ("Control" if t.endswith("Control") else "NA"))
    return M, probes, np.array(groups), np.array(samples)


def parse_annot(raw):
    txt = gzip.decompress(raw).decode("utf-8", "ignore")
    lines = txt.splitlines()
    begin = None
    for i, l in enumerate(lines):
        if l.startswith("!platform_table_begin"):
            begin = i
            break
    hdr = lines[begin + 1].split("\t")
    low = [h.lower() for h in hdr]
    id_idx = low.index("id")
    sym_idx = next(i for i, h in enumerate(low) if "gene symbol" in h)
    probe2sym = {}
    for l in lines[begin + 2:]:
        if l.startswith("!platform_table_end"):
            break
        p = l.split("\t")
        if len(p) <= max(id_idx, sym_idx):
            continue
        pid = p[id_idx].strip().strip('"')
        sym = p[sym_idx].strip().strip('"')
        if not sym or sym in ("?", "---", "NA"):
            continue
        sym = sym.split("//")[0].split("|")[0].strip().upper()
        if sym and re.match(r"^[A-Z0-9\-]+$", sym):
            probe2sym[pid] = sym
    return probe2sym


def build_genesets(tfs=None):
    gsets = {}
    with open(DOROTHEA, encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            src, tgt = parts[0].strip().upper(), parts[1].strip().upper()
            if src in (tfs or TFS):
                gsets.setdefault(src, set()).add(tgt)
    return gsets


def to_gene_matrix(M, probes, probe2sym):
    # collapse probes -> gene (mean across probes)
    from collections import defaultdict
    bygene = defaultdict(list)
    for j, pid in enumerate(probes):
        s = probe2sym.get(pid)
        if s:
            bygene[s].append(j)
    genes = sorted(bygene.keys())
    G = M.shape[1]
    expr = np.empty((len(genes), G), dtype=float)
    for i, g in enumerate(genes):
        idx = bygene[g]
        col = M[idx, :]
        expr[i] = np.nanmean(col, axis=0)
    # fill remaining nan with column median
    for j in range(G):
        c = expr[:, j]
        m = np.nanmedian(c)
        c[np.isnan(c)] = m
    return expr, genes


def ssgesva(expr, genes, gsets):
    """预计算 rank 矩阵 (genes x samples, ascending: 高表达->高 rank)，
    然后用 rank-mass fraction 给任意基因集打分: score = mean(rank in set)/G, 范围[0,1]。"""
    G = len(genes)
    N = expr.shape[1]
    R = np.empty((G, N))
    for j in range(N):
        R[:, j] = stats.rankdata(expr[:, j])
    gene_idx = {g: i for i, g in enumerate(genes)}

    def score_set(idx):
        return R[idx, :].sum(axis=0) / (len(idx) * G)

    out = {}
    for name, st in gsets.items():
        idx = np.array([gene_idx[s] for s in st if s in gene_idx])
        present = len(idx)
        out[name] = (score_set(idx), present, len(st))
    return out, R, G, N


def main():
    # 1. load / download bulk
    if os.path.exists(BULK_CACHE):
        d = np.load(BULK_CACHE, allow_pickle=True)
        expr, genes, groups = d["expr"], list(d["genes"]), list(d["groups"])
        print("[load] cached bulk:", expr.shape, "groups", dict(__import__("collections").Counter(groups)))
    else:
        print("[download] matrix + annot ...")
        M, probes, groups, samples = parse_matrix(download(MATRIX_URL))
        ann = parse_annot(download(ANNOT_URL))
        print("[annot] probe->symbol mapped:", len(ann))
        expr, genes = to_gene_matrix(M, probes, ann)
        np.savez(BULK_CACHE, expr=expr, genes=np.array(genes, dtype=object), groups=np.array(groups))
        print("[save] bulk cache:", expr.shape)

    print("[genes] total gene-space:", len(genes))
    # 2. gene sets
    gsets = build_genesets(TFS + NC_TFS)
    for t in TFS + NC_TFS:
        if t in gsets:
            print(f"  DoRothEA {t}: {len(gsets[t])} targets")

    # 3. ssGSEA (returns scores + rank matrix)
    scores, R, G, N = ssgesva(expr, genes, gsets)

    # 4. stats (with negative controls)
    is_stroke = np.array([g == "Stroke" for g in groups])
    is_ctrl = np.array([g == "Control" for g in groups])
    nS, nC = int(is_stroke.sum()), int(is_ctrl.sum())

    def score_idx(idx):
        return R[idx, :].sum(axis=0) / (len(idx) * G)

    def diff_of(sc):
        return float(np.mean(sc[is_stroke]) - np.mean(sc[is_ctrl]))

    rng = np.random.default_rng(2026)
    rows = []
    ps = []
    all_names = [t for t in (TFS + NC_TFS) if t in scores]
    for name in all_names:
        sc, present, total = scores[name]
        s, c = sc[is_stroke], sc[is_ctrl]
        U, p = stats.mannwhitneyu(s, c, alternative="two-sided")
        cliff = 2 * U / (nS * nC) - 1
        obs_diff = diff_of(sc)
        # random negative-control sets of matched size
        n_rand = 500
        rdiffs = []
        for _ in range(n_rand):
            ridx = rng.choice(G, size=present, replace=False)
            rdiffs.append(diff_of(score_idx(ridx)))
        rdiffs = np.array(rdiffs)
        emp_p = (np.sum(np.abs(rdiffs) >= abs(obs_diff)) + 1) / (n_rand + 1)
        rows.append({
            "module": name, "negative_control": name in NC_TFS,
            "n_targets_total": total, "n_targets_present": present,
            "stroke_mean": float(np.mean(s)), "control_mean": float(np.mean(c)),
            "mean_diff": obs_diff, "cliffs_delta": float(cliff),
            "U": float(U), "p": float(p), "empirical_p_random": float(emp_p),
        })
        ps.append(p)

    # BH q across the primary (non-control) modules
    prim = [r for r in rows if not r["negative_control"]]
    pp = np.array([r["p"] for r in prim])
    order = np.argsort(pp)
    q = np.empty_like(pp)
    m = len(pp)
    cummin = 1.0
    for rank, i in enumerate(order):
        cummin = min(cummin, pp[i] * m / (rank + 1))
        q[i] = cummin
    for r, qq in zip(prim, q):
        r["q_bh"] = float(min(qq, 1.0))

    print("\n=== ssGSEA module activity: Stroke vs Control (GSE16561, blood) ===")
    for r in rows:
        tag = "  [NEG-CTRL]" if r["negative_control"] else ""
        qstr = (f" q_BH={r['q_bh']:.2e}" if not r["negative_control"] else "")
        print(f"  {r['module']:8s} present={r['n_targets_present']:4d}/{r['n_targets_total']:4d}  "
              f"stroke={r['stroke_mean']:.4f} control={r['control_mean']:.4f}  "
              f"diff={r['mean_diff']:+.4f}  CliffΔ={r['cliffs_delta']:+.3f}  "
              f"p={r['p']:.2e}  empPRand={r['empirical_p_random']:.3f}{qstr}{tag}")

    json.dump({"dataset": "GSE16561", "tissue": "peripheral whole blood",
               "platform": "GPL6883 (Illumina HumanWG-6)",
               "n_stroke": nS, "n_control": nC,
               "method": "AUCell/ssGSEA-style rank-mass fraction activity score",
               "modules": rows},
              open(os.path.join(HERE, "human_bulk_gsva.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n[done] -> human_bulk_gsva.json")


if __name__ == "__main__":
    main()
