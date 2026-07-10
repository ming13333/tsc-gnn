# L5c 阶段摘要（Stage Summary · 2026-07-09）

> 本文件由 lean-context 技能压缩生成，供后续轮次复用，避免重复全量重读 v3 草稿/脚本。

## Session Intent
完成 L5 证据层的 **A+C 模式** 中 C 环节 = **公开 sc-CRISPR / perturb-seq 重分析**（Replogle 2022, K562 全基因组 CRISPRi pseudo-bulk），为"因果图从未被扰动"质疑补单细胞分辨率的方向性因果支持。
- A 环节（已完成）：公开 bulk TF-KO（GSE269122 Sox10 少突 cKO；GSE273163 Cebpb 杂合 KO 库普弗）
- B/L5b（已完成）：SigCom LINCS 基因扰动方向一致性（GATA2 OE rank 3/33782）
- C 环节（本次）：Replogle 2022 sc-CRISPR

## Files Created
- `gate1/rewiring_study/l5_perturbation/l5c_replogle.py` — L5c 主脚本：anndata backed 读 h5ad → 识别扰动基因列 → 提取目标 TF 行 + NTC 行 → 算 CRISPRi vs NTC 基线校正效应 → MWU(less)+Fisher OR(bottom-quartile)+rank among all candidate programs
- `gate1/rewiring_study/l5_perturbation/K562_gwps_normalized_bulk_01.h5ad` — Replogle 2022 K562 genome-scale pseudo-bulk（gemgroup Z-normalized，每扰动一行），357MB，Figshare 文件 ID 35773217，下载于 2026-07-09 完成（358MB）

## Files Read / Reused (no change)
- `gate1/data/dorothea/human_dorothea_regulon.tsv` — 人类 DoRothEA 靶集，格式 `source\ttarget\tweight\tconfidence`；SOX10=322 / CEBPB=589 / GATA2=5370 个靶
- `gate1/rewiring_study/Methods_Results_初稿_v3_2026-07-09.md` — v3 草稿（45KB），L5c 整合目标
- `gate1/rewiring_study/l5_perturbation/L5b_SigCom_报告_2026-07-09.md` — L5b 报告（风格参照）

## Decisions Made
- 数据源选 Replogle 2022 K562_gwps pseudo-bulk（357MB，每扰动一行）而非 65GB 单细胞矩阵或 48GB MTX → 避免加载百万级细胞
- 预期 **SOX10 在 K562（髓系白血病）不表达 → 该项诚实披露为 ABSENT**；CEBPB/GATA2 应存在扰动
- 位置：L5c 定位为单细胞分辨率 **causal support（程序级方向一致性）**，非 validation（K562 癌细胞系非脑、非卒中、CRISPRi 敲低非敲除）
- 红线：靶程序取自 DoRothEA GRN（非循环）；只验 TF/模块级方向，禁"validated/proves"

## Current State
- 下载：COMPLETE（358MB）
- 脚本：已写，待运行
- 环境：bbb_gnn 已装 anndata 0.11.4 / h5py 3.16.0 / mygene（var index 若为 Ensembl 则自动映射）
- 任务 #60 完成（可行性/TF覆盖）、#61 进行中（下载+提取）

## Next Steps
1. 运行 `l5c_replogle.py`（bbb_gnn env）→ 产出 `l5c_replogle_result.json`
2. 分析 CEBPB/GATA2 CRISPRi 程序级方向性富集 + TF×TF 交叉特异性；检查 SOX10 是否 ABSENT
3. 写 `L5c_Replogle_scCRISPR_报告_2026-07-09.md`（风格仿 L5b）
4. 整合进 v3 草稿：§2.5 证据阶梯表补 L5c 行；§3.7 新增 sc-CRISPR 子节；§4.2④ 与 §4.3 补单细胞分辨率诚实边界（SOX10 缺失/K562 癌细胞系/非卒中）；摘要；待办清单补 L5c check
5. 更新 `2026-07-09.md` 当日日志 + MEMORY.md（L5c 结果）

## Open Questions / Risks
- Replogle h5ad `obs` 扰动基因列名未知（脚本按 ["gene","perturbation","target_gene","guide_target","symbol","target","guide_id"] 探测）
- NTC 检测用行首锚定正则 `^(NON|NTC?|CONTROL|INTERGENIC|SAFE|NEG|EMPTY|MOCK|UNTREATED|UNKNOWN)`（避免 CNTN1 误判）；若无 NTC 则 fallback 到非 DoRothEA 基因
- var index 若为 Ensembl（ENSG 开头）则 mygene 映射；若为基因符号则直接用
- backed 模式只读目标行，避免全矩阵加载
