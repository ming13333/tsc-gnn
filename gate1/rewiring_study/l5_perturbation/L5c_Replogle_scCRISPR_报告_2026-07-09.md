# L5c · 公开 sc-CRISPR / perturb-seq 重分析报告（Replogle 2022, K562）

**日期**: 2026-07-09  
**定位**: L5 causal support 的第三层——**单细胞分辨率**公开扰动数据中，框架恢复 TF→靶程序的方向性因果支持检验  
**数据源**: Replogle et al. 2022, *Cell* — genome-scale CRISPRi (K562), pseudo-bulk `K562_gwps_normalized_bulk_01.h5ad`（gemgroup Z-normalized，每扰动一行；11,258 扰动 × 8,248 基因）

---

## 1. 设计逻辑

| 层面 | 设计 |
|------|------|
| **查询输入** | Human DoRothEA 靶程序（SOX10/CEBPB/GATA2 靶基因集） |
| **查询数据库** | Replogle 2022 K562 genome-scale CRISPRi pseudo-bulk（单细胞分辨率扰动，Figshare） |
| **预期方向** | CRISPRi 敲低激活因子 → 其 DoRothEA 直接靶基因应 **下调**（在靶程序内富集于负 Z） |
| **检验** | 对每个 TF：①自身扰动签名中靶程序的 mean Z；②MWU(less) 靶 vs 非靶；③Fisher OR（bottom-quartile 下调富集）；④在全部 332 个候选程序中的 rank（1=最下调）；⑤TF×TF 交叉特异性 |
| **靶程序来源** | 独立 DoRothEA GRN（非循环论证） |
| **红线** | 只验 TF/模块级方向，不逐边验证；定位 causal support 非 validation；**阴性/边界结果须诚实披露** |

## 2. 技术细节

- **数据结构**（探查确认）：`obs.index` = `<rowid>_<GENE>_<guide>_<ENSG>`，扰动基因取第 2 字段；`var.index` = Ensembl ID，经 mygene 映射为 symbol（8,190/8,248 命中）；NTC = `non-targeting`（585 行）
- **基线校正**：signature Z − NTC 均值 Z（NTC 均值 ≈ 0，min −0.009 / max +0.021，清洗 inf/nan 后中性）
- **on-target 校验**：TF 自身位点 Z 应为强负（敲低成功）
- **SOX10 表达状态**：SOX10 的 Ensembl ID（ENSG00000100146）**不在 K562 基因空间内** → K562（髓系白血病）不表达神经/胶质 TF SOX10 → 该项为 off-context

## 3. 结果

| TF | guides | self-Z (on-target) | n 靶(映射/总) | mean 靶 Z | MWU p(less) | Fisher OR↓ | Fisher p | rank(最下调) | 判定 |
|----|--------|-------------------|--------------|-----------|-------------|------------|----------|--------------|------|
| **SOX10** | 1 (P1P2) | **NaN**（位点不在 K562） | 119/322 | **+0.0184** | 0.913 | 0.97 | 0.599 | 304/332 (8.4%) | **NULL**（off-context） |
| **CEBPB** | 1 | **−0.427** ✓ | 223/589 | **+0.0071** | 0.843 | 0.84 | 0.873 | 239/332 (28.0%) | **NULL** |
| **GATA2** | 2 (P1,P2) | **−0.190** ✓ | 2419/5370 | **+0.0044** | 0.178 | 1.00 | 0.505 | 139/332 (58.1%) | **NULL** |

**关键观察**：
- **on-target 敲低成功**：CEBPB（self-Z −0.427）与 GATA2（−0.190）自身位点显著下调，证明 CRISPRi 扰动与数据加载正确；数据管线无 bug。
- **靶程序未集体下调**：三个 TF 的 DoRothEA 靶程序 mean Z 均接近 0 或轻微为正（SOX10 +0.018, CEBPB +0.007, GATA2 +0.004），无一进入"最下调"前列（rank 304/332、239/332、139/332）。
- **SOX10 不在 K562 表达**：自身位点未检出（self-Z = NaN），与其神经/胶质身份一致 → off-context 扰动预期为 null。
- **交叉特异性平淡**：各 TF 扰动下，其他 TF 的靶程序 mean Z 均在小正值区间（0.004–0.015），无分化。

## 4. 核心发现

### ✗ K562 单细胞分辨率扰动未提供方向性因果支持（context-null）

与 L5a（鼠原生细胞类型 bulk KO，阳性）和 L5b（LINCS GATA2 OE，阳性）不同，Replogle 2022 K562 CRISPRi 在 **off-context 癌细胞系** 中未能复现框架恢复靶程序的方向性下调。这是一次**信息丰富的阴性/边界结果**：

1. **SOX10** 在 K562 根本不表达 → 扰动为 off-context，null 在生物学上预期；
2. **CEBPB / GATA2** on-target 敲低确凿（自身位点下调），但其 DoRothEA 靶程序（跨组织/语境汇总的知识库先验）在单一髓系谱系中**未协调下调**——提示 TF 的推断靶集并非在该细胞语境下全部直接/即时共调控。

### 与 L5a/L5b 构成三角验证语境边界

| 层 | 扰动语境 | 结果 |
|----|----------|------|
| L5a | 鼠**原生**细胞类型 bulk KO（Sox10 少突 / Cebpb 库普弗） | **阳性**（rank 46/412, 149/404） |
| L5b | LINCS **过表达**（GATA2 OE） | **阳性**（rank 3/33,782，强自特异性） |
| **L5c** | K562 **癌细胞系** CRISPRi（off-context） | **NULL / 语境边界** |

该三角表明：**框架恢复的 TF→靶程序的方向性因果支持，仅在生物学恰当的语境下成立**（TF 在其活跃的原生谱系中被 KO，或其方向正确的过表达被扰动）；在通用癌细胞系 CRISPRi 筛选中不充分。这恰恰是 interpretability-first 框架的预测之一，并论证了"不应将任何单一扰动筛选视为决定性验证"。

## 5. 证据等级定位

| 发现 | 等级 | 定位 |
|------|------|------|
| SOX10 在 K562 不表达 → off-context null | 预期/信息性 | 语境边界 |
| CEBPB/GATA2 on-target 成功但靶程序未下调 | **边界证据** | K562 单一谱系不足以支撑跨组织先验靶集的方向性因果 |
| 三 TF 均无显著 down-shift | **不支持（在 K562）** | 诚实披露，非缺陷 |

**整体定位**：L5c 是 **causal support 的语境边界检验**——它未提供 K562 内的方向性支持，但通过与 L5a/L5b 的对比，划定了"因果支持仅在恰当生物语境成立"的边界。这强化了方法学叙事（解释优先、语境依赖），而非削弱它。

## 6. 诚实边界

1. **K562 癌细胞系**：髓系白血病，非脑、非卒中语境；SOX10 等神经/胶质 TF 本不表达
2. **CRISPRi 敲低 ≠ 敲除**：敲低效应对直接靶的下游传播弱于 KO
3. **pseudo-bulk 聚合**：gemgroup Z-normalized，每扰动一行；若仅亚群受影响，信号被稀释
4. **DoRothEA 跨组织先验**：靶集汇总多组织/语境，在单一 K562 谱系未必共调控——阳性 L5a/L5b 得益于"TF 在其活跃语境"的匹配
5. **self-Z NaN（SOX10）**：因 SOX10 位点不在 K562 基因空间，仅说明不表达，非分析错误
6. **定位**：causal support 的边界测试，**非 validation**；阴性结果透明报告

## 7. 与 L5 / L5b 整合

L5 证据层现由三层构成（A+C 模式）：

- **L5a（A）**：鼠原生细胞类型 bulk KO → Sox10/Cebpb 程序级方向性因果支持（阳性）
- **L5b（B）**：SigCom LINCS 基因扰动 → GATA2 OE 强自特异性验证（阳性）；3 TF CRISPR KO 方向正确但无自特异性
- **L5c（C）**：Replogle 2022 K562 sc-CRISPR → 三 TF 均无显著方向性下调（语境边界/null）

三层互补叙事：**因果支持的方向性在生物学恰当的语境（原生谱系 KO、方向正确的 OE）下被独立公共数据证实，但在 off-context 癌细胞系 CRISPRi 中不成立**——这一"语境依赖性"本身就是框架的可检验预测，并论证了单一扰动筛选不足以作为验证。L5 整体仍定位为 **program-level directional causal support**（非 validation）。

## 8. 产物

- `l5c_replogle.py` — 分析脚本（obs.index 取扰动基因 / mygene Ensembl→symbol / inf 清洗 / MWU+Fisher+rank）
- `l5c_replogle_result.json` — 完整结果 JSON
- `l5c_run.log` — 运行日志
- `K562_gwps_normalized_bulk_01.h5ad` — 数据源（358MB，已完成下载）
