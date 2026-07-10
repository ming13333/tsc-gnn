# 虚拟敲除项目 · 会话主压缩摘要（Session Master Summary）

> lean-context 压缩产物（2026-07-09，cycle 2）。后续轮次优先复用本文件 + `MEMORY.md`，**勿重复全量重读 v3 草稿(45KB)/脚本/原始日志**。
> 关键事实以 `MEMORY.md` 为权威；本文件补"会话级产物轨迹 + verbatim 标识符 + 数值"。

## Session Intent
构建纯计算方法学框架"虚拟扰动(Virtual Perturbation)"，冲 Patterns / Genome Biology / Nat Methods / Bioinformatics，无湿实验。
主轴（已定调，勿重议）：图不优于线性预测 → 卖点 = **边级图重布线可解释性(ΔW_ij + 置换显著性)** + 时间×状态×药物逆转。对标 CellOracle/SCENIC/NicheNet。
当前阶段：L1–L5 五层证据已全部落地，A+C 模式（L5 公开扰动因果支持）已闭合。下一步 = 全文 polish / 选刊。

## 关键决策（verbatim）
- 不做应用文，做**方法框架文**。
- 深核 = TSC-GNN（时间×状态条件化 GNN，GRN+通讯图重排）；"TSC-GNN 超越线性"被真实验证证伪（卒中 held-out −22.3%；GEARS 单 −8.3%/双 −4.2%，CI 含 0）→ 严禁作主卖点。
- 鲁棒性强结论：10 种子 × 5 图(k-NN/DoRothEA/random/perm-DorothEA/0-hop) × 2 任务 = 90 点；图优于线性 **0 例**、无差异 86 例 → "固定图+线性读头图未带来预测优势"为强结论。
- L5 红线：= **causal support 非 validation**（重分析公开数据 / 非卒中语境 / 模块级非边级 / 非自有扰动）；阴性须诚实披露为语境边界，非缺陷。

## 证据阶梯 L1–L5（落地于 v3 草稿）
- **L1 技术复现(鼠内跨队列)**：TF 主调控因子排序 ρ=0.48–0.55 (p<1e-15)；Sox10 三队列 Top20 共通；边级不重现(Jaccard≈0) → 可解释性主轴卖点。
- **L2 已知生物学重现**：Sox10 等已知 TF 在 target 富集（OR 36–61 少突/神经元 OR 110/小胶 OR 15；14/15 已知 TF 出现）。
- **L3 跨物种调控汇聚**：C2 人 GRN 结构 SOX10→髓鞘 OR=27 emp_p=0.0005；CEBPB→炎症 OR=16.5 p=0.024；GATA2→炎症 OR=4.6 p=0.047。+ GSE16561 人卒中 bulk 激活 AUCell/ssGSEA BH-q=1.1e-8（PAX6 负对照也激活→非唯一）。披露外周血非脑/横断面/未反卷积。
- **L4 药物逆转 LINCS**：L1000CDS2(maayanlab JSON 免 key)。robust 命中 vorinostat/trichostatin A(HDAC 抑制剂)+Mevastatin/Rosuvastatin(他汀)；置换 real=4 vs 随机 2.75±1.37 max5, p=0.33 **不显著** → hypothesis generation 非药效预测。vorinostat 跨系一致是 L1000 背景（53% 随机签名出现）。
- **L5 公开 TF 扰动方向性因果支持（A+C 模式，已全部完成）**：
  - **A（bulk KO）**：GSE269122 鼠少突 Sox10 cKO → Sox10 靶程序 rank 46/412 前11% 最下调 OR↓=1.81 p≈0；GSE273163 鼠 Cebpb 杂合 KO 库普弗 → Cebpb 靶程序 rank 149/404 前37% OR↓=1.20 p=0.022，且 Sox10 靶在此不下调(OR=0.95 p=0.68) = 干净特异性。
  - **B / L5b（SigCom LINCS）**：GATA2 OE rank 3/33782 top0.01% mimicker (p=1.4e-05，强自特异性 self 0.01% vs cross 4.98%)；3 TF CRISPR KO 全 top~1% reverser（方向对但无 TF 级自特异性）；SOX10 OE 方向混合；CEBPB OE 无数据。
  - **C（sc-CRISPR / L5c，本次完成）**：Replogle 2022 K562 genome-scale CRISPRi pseudo-bulk。三 TF 靶程序在自身 CRISPRi 下**均不显著下调** → 语境边界 null（详见下"Files Created / L5c 结果"）。与 A（原生谱系 KO 阳性）+ B（GATA2 OE 阳性）构成**三角验证**：「因果支持仅在生物学恰当语境成立」——此语境依赖恰为 interpretability 框架的可检验预测，强化 causal support 非 validation 定位。

## Files Created（本会话产物，verbatim 路径）
- `gate1/rewiring_study/Methods_Results_初稿_v3_2026-07-09.md` — v3 草稿(45KB)，方法+结果全文，L5c 已整合(§2.5/§3.7/§4.2④/§4.3/摘要/待办/源产物)。
- `gate1/rewiring_study/l5_perturbation/l5c_replogle.py` — L5c 主脚本（anndata `backed='r'` 读 h5ad；扰动基因取 `obs.index` 第 2 个 `_` 字段 `split('_')[1].upper()`；NTC 精确匹配 `'non-targeting'`；`mannwhitneyu(alternative='less')` + `fisher_exact(alternative='greater')`；rank = `(all_means < self_mean).sum()+1`；mygene 映射 Ensembl var；`np.nan_to_num` 清 inf/nan）。
- `gate1/rewiring_study/l5_perturbation/K562_gwps_normalized_bulk_01.h5ad` — Replogle 2022 K562 pseudo-bulk（gemgroup Z，每扰动一行），358MB，Figshare 文件 ID **35773217**（article 20029387）。
- `gate1/rewiring_study/l5_perturbation/l5c_replogle_result.json` — L5c 数值结果。
- `gate1/rewiring_study/l5_perturbation/L5c_Replogle_scCRISPR_报告_2026-07-09.md` — L5c 报告（语境边界 null + 三角验证）。
- `gate1/rewiring_study/l5_perturbation/L5c_stage_summary.md` — ⚠️ 陈旧（写于分析运行前，"脚本待运行"），以本 Master Summary 为准。
- 历史产物（已存在，未本会话改）：`L5_upgrade_公开扰动因果支持_2026-07-09.md`、`L5b_SigCom_报告_2026-07-09.md`、`CROSS_COHORT_报告_2026-07-09.md`、`CROSS_SPECIES_报告_2026-07-09.md`、`L3_upgrade_人bulk激活保守_2026-07-09.md`、`StepD_文章级分析阐述_2026-07-09.md`、`Validation_章草稿_2026-07-09.md`、`AUDIT_全流程审查_2026-07-09.md`。

## L5c 实测数值（verbatim，来自 l5c_replogle_result.json）
- 数据规模：11,258 扰动 / 585 NTC（`non-targeting`）。
- **SOX10**：自身位点 ENSG00000100146 **不在 K562 var 基因空间**（髓系白血病不表达神经/胶质 TF）→ off-context；靶程序 mean Z=+0.018，rank 304/332（非下调）。self-Z=NaN。
- **CEBPB**：self-Z=−0.427（on-target 敲低成功）；靶程序 mean Z=+0.007，MWU p=0.84，Fisher OR 不显著，rank 239/332。
- **GATA2**：self-Z=−0.190（on-target 成功）；靶程序 mean Z=+0.004，MWU p=0.18，rank 139/332。
- 解释：K562（癌细胞系）非神经/胶质语境，CEBPB/GATA2 虽 on-target 但靶程序未在单细胞分辨率集体下调；与 A（原生谱系 bulk KO 阳性）形成语境边界对比。

## Files Modified（本会话）
- `gate1/rewiring_study/Methods_Results_初稿_v3_2026-07-09.md` — 整合 L5c：§2.5 证据阶梯表补 sc-CRISPR 行；§3.7 新增 Replogle 子节；§4.2④ 补三角验证；§4.3 补语境边界诚实披露；摘要补语境依赖；待办清单新增 L5c 完成项；源产物清单加入 L5c 文件。
- `.workbuddy/memory/MEMORY.md` — L5 段补 C(已完成,语境边界 null) + 三角验证要点（line 32–35）。
- `.workbuddy/memory/2026-07-09.md` — 当日日志追加 L5c 完成记录。

## 运行环境（非显然，必记）
- 解释器 = conda `bbb_gnn`：`C:/miniconda3/envs/bbb_gnn/python.exe`（numpy 2.2.6 / scipy 1.15.3 / pandas 2.3.3）。managed venv 与 base conda 均无 numpy → **必须用 bbb_gnn**。
- 已装：`anndata 0.11.4` / `h5py 3.16.0` / `mygene`（L5c 用）。
- 长任务前台 + `tee` + `timeout 600000`；后台 shell 会话挂起会被回收。

## Current State
- L1–L5 全部落地；A+C 模式 L5 证据链闭合（A 阳性 / B 阳性 / C 语境边界 null）。
- v3 草稿已含 L5c 整合；MEMORY.md / 当日日志已更新。
- 任务 #60/#61/#62/#63 全部 completed。

## Next Steps（建议，待用户确认）
1. 全文 polish：abstract 重读、图注统一、缩写表、参考文献格式。
2. 选刊评估：Patterns（方法学友好）/ Genome Biology / Nat Methods / Bioinformatics / Cell Systems；比对 ARTEMIS2025/RegVelo2026 等对手定位差异化。
3. 补充图形：L5 三角验证示意图、时间轴×状态×药物逆转闭环图。
4. 预印本（bioRxiv）或直投（依选刊）。

## Open Questions / Risks
- L4 置换不显著(p=0.33)须诚实定位 hypothesis generation。
- L3 外周血非脑/横断面/PAX6 负对照须披露。
- 基础模型威胁（scGPT/scFoundation/Geneformer）需在 Discussion 论证差异化（时间×状态+可解释性）。
- DoRothEA 为因果方向图但源于 bulk，rewiring 解释可信度前提已建。
- 跨物种 HRA007397 实为受控访问 PBMC fastq（非脑）→ 已改正交投影方案。
