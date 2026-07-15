# 运算过程符合性审计（Conformance Check）— manuscript_v7

**日期**：2026-07-14 ｜ **方法**：repro-guard 独立回读脚本/CSV/JSON，不采信记忆；严格对齐口径分母
**范围**：manuscript_v7.docx（= manuscript_v7.md 同步版）全部数值结论 vs 底层计算产物

---

## 0. 总判定

| 维度 | 状态 | 说明 |
|---|---|---|
| **数值结论一致性** | 🟢 **全部对齐** | 稿件所有关键数字（L1–L5、3.3、GEARS、90-config、L4）均与底层 CSV/JSON/报告逐一对齐，无数值矛盾 |
| **产物/呈现一致性** | 🔴 **1 处严重 + 1 处中等** | Fig2A 合成散点但图注称"真实 90 配置散点"；L5b 结果 JSON 被失败重跑覆盖为空 |

**一句话**：数字本身站得住（比 0→4 时期扎实得多）；问题在**两张图/产物的可复现性与披露**，而非算错。最严重的是 Fig2A 的图注与实现不符。

---

## 一、逐项核验（✅ = 稿件数字 = 计算产物）

| 稿件结论 | 来源产物 | 核验值 | 判定 |
|---|---|---|---|
| 90-config：4 better / 0 worse / 86 n.s.（全 seed-1 combo） | `gears_validation/robustness_results.csv` | 逐行数 `sig` 列：ridge 图变体（排除 +deg/kernel）恰好 4 better（seed1 combo: knn/doro/rand/0-hop）+ 0 worse + 86 n.s. | ✅ |
| GEARS：−8.3%（single）/ −4.2%（double），CI 含 0 | `gears_validation/run_gears.py`、`RESULTS.md` | rel_imp=−8.3% CI[−76%,+34%]；−4.2% CI[−27%,+13%] | ✅ |
| 3.3 边计数 8/28/19/36；ΔW +0.51/+0.59/+0.79/+0.42 | `REWIRING_报告_v3_2026-07-09.md` | pooled-q<0.05：8/28/19/36；Sox10→Plp1 2d→14d +0.51、Hey2→Acta2 +0.59、Sox9→Hapln1 +0.79、Sox10→Plp1 sham→14d +0.42 | ✅ |
| L1：ρ=0.48–0.55，p<1e-15，Jaccard≈0，dir≈0.52 | `CROSS_COHORT_报告_2026-07-09.md` | ρ=0.517/0.548/0.482，p=1.5e-18/8.8e-20/1.7e-15；方向一致≈0.52；Jaccard≈0 | ✅ |
| L3a：OR 27/16.5/4.6，emp_p 0.0005/0.024/0.047 | `human_module_enrich.json` | SOX10→myelin OR=27.00 emp_p=0.00050；CEBPB→neuro OR=16.49 emp_p=0.0235；GATA2→neuro OR=4.62 emp_p=0.0465 | ✅ |
| L3b：BH-q=1.1e-8，emp_p=0.002，PAX6 p=4.5e-8 | `human_bulk_gsva.json` | q_bh=1.12e-8；empirical_p_random=0.0020；PAX6 p=4.49e-8 | ✅ |
| L4：0.0671 / 0.125 / 4 hits / 2.75±1.37 / p≈0.33 | `drug_reversal_result.json` + `drug_perm_result.json` | main best 0.0671；robust best 0.125；real_known_hits=4；perm_mean=2.75 std=1.37 max=5 p=0.333 | ✅ |
| L5a：rank 46/412；149/404，OR 1.20/0.95，p 0.022/0.68 | `l5_perturbation/l5_causal_direction.json` | GSE269122 Sox10 rank46/412；GSE273163 Cebpb rank149/404 OR=1.198 p=0.022；Sox10 OR=0.953 p=0.682 | ✅ |
| L5b：GATA2 OE rank 3/33782，p=1.4e-5，self 0.01% vs cross 4.98% | `L5b_SigCom_报告` + `l5b_run.log`（首次成功段） | GATA2→GATA2 rank3 mimicker top0.01% p=1.4e-5；best cross (GATA2→SOX10) rank32101=top4.98% | ✅* |
| L5c：self-z −0.43/−0.19，meanZ +0.007/+0.004，MWU 0.84/0.18，rank 239/139 | `l5_perturbation/l5c_replogle_result.json` | CEBPB self_z=−0.427 meanZ=0.00706 mwu=0.843 rank239；GATA2 self_z=−0.190 meanZ=0.00437 mwu=0.178 rank139 | ✅ |
| 3.8：MYC 19/332 p=3.1e-3；BCL11A 29/332 p=7.5e-3；GATA1 rank187 | `l5_perturbation/l5c_positive_control_result.json` | MYC rank19 mwu=0.00314；BCL11A rank29 mwu=0.00754；GATA1 rank187 meanZ=0.268 | ✅ |

\* L5b 数字对，但**结果 JSON 是空的**（见下文 B1）。

---

## 二、不一致 / 风险环节（按严重度）

### B1. 🔴【严重】Fig2A 是合成散点，图注却称"真实 90 配置散点" —— 披露/诚信问题

- **图注（manuscript_v7.md L372）**："(A) Scatter plot of graph-vs-linear relative improvement across 90 configurations."
- **实现（make_figures.py L14–17, L81, L95）**：
  - L14–17：*"Fig 2A individual points are a transparent reconstruction of the audited 4/86/0 outcome (**raw per-configuration values were not archived in this snapshot**)"*
  - L81：`rng = np.random.default_rng(20260710)`
  - L95：`y = rng.normal(0.0, 0.9)` —— 无差异云是**随机合成**的；4 个 better 点是硬编码偏移。
- **矛盾**：图注断言画的是 90 个真实配置，实际是随机点 + 4 个手工点。文字结论（4/86/0）本身由 `robustness_results.csv` 支撑（见上表 ✅），但**图没有描绘真实数据**。
- **风险**：reviewer 要求原始散点即穿帮；属于"图注与实现不符"的诚信披露问题，非数字错误。
- **修法（推荐）**：`robustness_results.csv` 每个配置都有真实 `rel_imp` 与 `ci_lo/ci_hi` 列，**直接用 CSV 重画 Fig2A**（真实散点 + 4 个 CI 不含 0 的 better 点），删掉 `rng.normal` 合成逻辑。若坚持保留示意，至少把图注改为 "schematic reconstruction of the 4/86/0 summary"。

### B2. 🟠【中等】L5b 结果 JSON 被失败重跑覆盖为空

- **证据（l5b_run.log）**：脚本首次成功算出 3×3 矩阵（L19–35，含 GATA2 OE rank 3），L37 写盘 `[done]`；随后第二次重跑 GATA2 Overexpression 查询 **HTTP 500 失败 8 次**（L40–56），矩阵全 `absent`，L61 再次 `[done]` **覆盖**了 JSON。
- **现状**：磁盘 `l5b_sigcom_result.json` 是失败重跑的空壳（矩阵全 `[]`，仅存 targets_n）。稿件数字来自首次成功（log L19–35 / `L5b_SigCom_报告`），**数字正确，但可机读产物坏掉**。
- **风险**：重跑脚本会复现失败→得到空 JSON；结果不可持久复现。
- **修法**：① 从成功 log 段重生成干净 JSON 落盘；② 脚本改为"仅全部查询成功才写盘"+ 对 500 做退避重试 + 不二次覆盖。

### B3. 🟡【轻微】GEARS 措辞 "−8.3% … improvement" 反直觉

- `rel_imp=−8.3%` 意为"图比线性差 8.3%"，并非 improvement。稿件 3.1 写 "gave −8.3 % … improvement, both with confidence intervals including zero"——数字对（−8.3/−4.2、CI 含 0），但 "improvement" 与负号冲突。
- **建议**：改 "no significant gain (graph −8.3% vs linear, CI including zero)"。

### B4. 🟡【轻微】kernel 探针计数口径

- CSV 实际含 **12 个 kernel 行全 worse**（seed 0/4/8 各 4 图，含 +deg）；审计与稿件称"6 行"。若稿件"6 行"指 ridge-kernel 子集（不含 +deg）则成立，否则多算一倍。不影响主基准 4/86/0。
- **建议**：注明口径（ridge-kernel 6 行 vs 含 +deg 12 行）。

### B5. 🟡【轻微】方向保留率 69–79% 来源待确认

- v7 写 "69–79 % of significant edges (pooled 74 %) preserved direction"；旧稿 `Methods_Results_初稿` 写 "77–83%"。两者在 v7 内部自洽（Table 1 用 OR 15–110 汇总），但数值从 77–83 变 69–79 应可追溯至当前数据集重算，避免静默改动。
- **建议**：确认 69–79 来自当前 `sanity_check.py` 输出并标注。

---

## 三、范围与说明

- **CellChat（§2.4）**：今日新生成 `cellchat_rewiring.csv` / `cellchat_rewiring_sig.csv` / `cellchat_py.py` 等，但 **manuscript_v7.docx 仍是 13:34 版，未含 CellChat 结果**（见前次逻辑闭环评审的硬伤）。故本轮"运算结论符合"范围不含 CellChat（v7 中无其数字）；待其加入稿件时须单独做此符合性检查（尤其配体-受体显著性阈值与 Fig  caption）。
- **GNN 引擎角色**：所有报告量均为 ΔW（相关性差分）+ 置换 + 模块富集，GNN 传播产出的 latent embedding 与 ΔW 的关系未在正文展示（前次评审已提）；本轮未新增证据，维持原判断。
- 全部 L1–L5 与 3.3 数字均经**独立产物交叉验证**，可靠，无需重算。

---

## 四、建议优先级

1. **【阻断·诚信】** Fig2A 改用 `robustness_results.csv` 真实 rel_imp 重画，或图注明确标注为 schematic reconstruction。（B1）
2. **【尽快·可复现】** L5b 重生成干净 JSON + 脚本加固防覆盖/防 API 失败。（B2）
3. **【顺手】** GEARS 措辞（B3）、kernel 计数口径（B4）、69–79% 来源确认（B5）。

---

## 五、修复记录（2026-07-14 晚）

### ✅ B1 已修复 — Fig2A 重画为真实数据
- `make_figures.py` 的 `fig2()` panel A 改为**直接读取** `../gears_validation/robustness_results.csv`：
  - 过滤条件 `model ∈ {gnn_knn, gnn_doro, gnn_rand, gnn_perm_doro, gnn_knn_nograph}` 且 `readout == ridge`（即 90-config 主基准：single 5 图 ×10 + combo 4 图 ×10 = 90）。
  - 断言 `len(rows) == 90`；逐点用真实 `rel_imp` 作 y，4 个 `better` 配置用真实 `ci_lo/ci_hi` 画非对称误差棒。
  - 删除了 `rng.normal(0,0.9)` 合成云与硬编码 better 点；抖动随机种子只用于水平散开。
  - 图注文字改为动态统计（"4 / 90 beat · 86 not significant · 0 worse"）。
- 核验：重跑脚本断言通过；90 点分布 4 better / 86 n.s. / 0 worse；rel_imp 范围 [-72.3, +19.1]；4 better 真实 CI 全部不含 0（如 seed1 combo k-NN +19.13 [+8.52,+29.60]）。
- 同步产物：`figures/figure2_prediction_benchmark.png`（位图）+ `.svg`（116 `<path>` / 28 `<text>`，可编辑矢量）。
- 已用 `build_docx_v7.py` 重建 `manuscript_v7.docx`，将真实 Fig2A 嵌入稿件；主副本与 `tsc-gnn-repo` 镜像均已同步。

### ✅ B2 已修复 — L5b JSON 重生 + 脚本加固
- 新增 `l5_perturbation/regen_l5b_json.py`：从 `l5b_run.log` 首次成功段（lines 1–37）**权威恢复**被失败重跑覆盖的真实结果（未再调 API，因彼时 SigCom 返回 500）。带断言：GATA2 OE→GATA2 rank 3 / 33782、p=1.4e-5、cross(GATA2→SOX10) rank 32101。
- 在两个副本（`gate1/rewiring_study` 与 `tsc-gnn-repo/gate1/rewiring_study`）均重生成 `l5b_sigcom_result.json`：18 个矩阵格中 13 个已填充（符合首次成功段的真实输出），含稿件所引 GATA2 OE rank 3/33782 头条。
- 加固 `l5b_sigcom.py`（两副本同步）：
  - 引入 `failed_any` / `n_failed` 标志；任一 query 失败（HTTP 500 等）时在**主输出路径之外**写 `.partial` 文件。
  - 改为**原子写**（临时文件 + `os.replace`），绝不留半截/损坏 JSON。
  - 关键：**失败重跑不再覆盖已有的好结果**——仅当全部 query 成功才写 `l5b_sigcom_result.json`，否则留警示并写 `.partial`。待 API 健康后重跑即得新鲜真值。
- 语法检查 `py_compile` 通过；两副本 JSON 经校验有效且头条数字正确。

### 剩余（非阻断）
- B3/B4/B5 仍为轻微措辞/口径项，未动；可视需要顺手修。
