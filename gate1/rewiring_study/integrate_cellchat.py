"""轻量整合: 读 cellchat_rewiring.csv + DoRothEA TSV,
计算"显著重布线的 LR 对"其配体/受体是否落于主调控因子(Sox10/Cebpb/Gata2)靶集,
并输出整合报告 + cellchat_rewiring_sig.csv。秒级运行, 不重算 MWU。"""
import os
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE
PROJ = os.path.dirname(os.path.dirname(HERE))
TSV = os.path.join(PROJ, "gate1", "data", "dorothea", "mouse_dorothea_regulon.tsv")
MRS = ["Sox10", "Cebpb", "Gata2"]


def load_targets():
    res = {}
    with open(TSV) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            res.setdefault(p[0], set()).add(p[1])
    return res


def main():
    targets = load_targets()
    mr_targets = set().union(*[targets.get(t, set()) for t in MRS])
    print(f"DoRothEA: {len(targets)} TF; 主调控因子靶集大小: "
          f"Sox10={len(targets.get('Sox10', set()))}, Cebpb={len(targets.get('Cebpb', set()))}, "
          f"Gata2={len(targets.get('Gata2', set()))}; 并集={len(mr_targets)}")

    df = pd.read_csv(os.path.join(OUT, "cellchat_rewiring.csv"))
    df["linked_to_MR"] = df.apply(
        lambda r: (r.ligand in mr_targets) or (r.receptor in mr_targets), axis=1)
    df.to_csv(os.path.join(OUT, "cellchat_rewiring_sig.csv"), index=False)

    sig = df[df.sig]
    print(f"\n显著 LR-转变记录: {len(sig)}")
    print(f"其中 配体/受体 落于 Sox10/Cebpb/Gata2 靶集: {int(sig.linked_to_MR.sum())} "
          f"({100*int(sig.linked_to_MR.sum())/max(1,len(sig)):.0f}%)")

    linked = (sig[sig.linked_to_MR][["ligand", "receptor", "pathway", "transition",
                                     "study", "log2FC", "sender", "receiver"]]
              .drop_duplicates(subset=["ligand", "receptor"]))
    print(f"\n=== 与主调控因子程序关联的 重布线 LR 对 (去重) ===")
    print(linked.to_string(index=False))

    # 跨转变一致的 top pairs (>=2 transitions 显著)
    pair_counts = sig.groupby(["ligand", "receptor"]).size()
    consistent = pair_counts[pair_counts >= 2].sort_values(ascending=False)
    print(f"\n=== 跨 >=2 个转变一致的 重布线 LR 对: {len(consistent)} ===")
    print(consistent.to_string())

    # 与 MR 关联的、且跨转变一致的 (最值得写进稿子的核心发现)
    cons_set = set(consistent.index)
    core = [p for p in cons_set if p in set(zip(linked.ligand, linked.receptor))]
    print(f"\n=== 核心: 既跨转变一致、又关联主调控因子程序的 LR 对: {len(core)} ===")
    for l, r in core:
        print(f"  {l} -> {r}  ({pathway_of(l, r)})")


def pathway_of(l, r):
    from cellchat_lrdb import pathway_of as pw
    return pw(l, r)


if __name__ == "__main__":
    main()
