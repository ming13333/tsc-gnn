# Gate 1 流水线审计报告（research-bug-audit · 模式 A）

- **审计日期**：2026-07-09
- **审计对象**：`gate1/` 全部代码与一次性实验结果（下载→整合→预处理→任务构造→评估）
- **审计目标**：校验既往流程**正确**（无静默 bug / leakage / 假阳性）且**可重复**（结果确定性、数据完整、环境可锁定）
- **结论**：流水线**逻辑正确、结果可重复、真实 FAIL 为真阴性（非假阴性）**。审计中发现 1 个阻断级隐患（已修复）+ 2 个重要级 + 3 个建议级，详见下表。

---

## 一、流水线全景（已逐一核对）

| 阶段 | 模块/函数 | 输入 | 输出 | 关键参数 |
|------|-----------|------|------|----------|
| 下载 | `data_acquisition.py` / `clean_download.py` / `gse225948_download.py` | GEO HTTPS | `data/GSE174574/*`(10x mtx), `data/GSE225948/*`(csv) | 已落盘，56/56 gzip 完整 |
| 整合 | `preprocessing.load_integrated_timeseries` | 两 GSE raw | 24h→2d→14d 轴 AnnData | cohorts={male,W*,MCAO} |
| 预处理 | `preprocessing.prep_for_gate` | 整合 adata | 2000 HVG + z-score + 状态 score | cell_ranger HVG, zscore 批次校正 |
| 任务 | `task_builder.build_heldout_task`(**方案A**) | prep 后 adata | held-out 基因表达预测任务 | holdout_frac=0.2, seed=0 |
| 评估 | `baselines` + `evaluate.evaluate_gate` | X/ctx/y/splits | rel_imp, bootstrap CI, corr, verdict | Ridge(α=1), 1000 bootstrap |

> 注：`build_real_task`(质心位移) 已被证明 leakage 并**弃用**（见阻断级 #1）。

---

## 二、可重复性实证（关键证据）

| 检查项 | 方法 | 结果 |
|--------|------|------|
| 合成信号 | 重跑 `generate_synthetic(seed=0)` | rel_imp=**94.6%** CI 93.7–95.4% PASS（与历史一致） |
| 合成反例 | 重跑 `generate_null(seed=1)` | rel_imp=**−1.0%** CI −1.2~−0.8% FAIL（与历史一致） |
| 真实 FAIL | 重跑 `run_gate1_timeseries.py` | MSE_lin=0.8819 / cond=0.8818 / rel_imp=**0.0%** CI 0.0–0.0% / corr=0.293 / **FAIL**（逐位一致） |
| 数据完整性 | `gzip.open` 读 56 个 .gz | **全部完整**，无截断/损坏 |
| 环境锁定 | 捕获版本 | numpy 2.4.6 / scipy 1.17.1 / pandas 3.0.3 / sklearn 1.9.0 / scanpy 1.12.2（已写入 `requirements_pin.txt`） |
| 正对照（敏感性） | 对 held-out 目标注入 state 驱动信号 | signal=0.5→3.2%; **signal=1.0→11.6% PASS** |

**正对照含义（最重要）**：真实 FAIL 不是任务「瞎」导致的假阴性——注入中等强度状态信号后条件化模型即显著恢复（rel_imp 11.6% > 10% 阈值）。因此真实 rel_imp≈0% 是**可信的真阴性**：在「共表达基线之外，炎症状态对 held-out 基因表达几乎无额外预测力」。

---

## 三、红 flag 静默 bug 扫描

| 红 flag | 结果 | 判定 |
|---------|------|------|
| 随机种子缺失 | `random_splits(seed=1)`、`build_heldout_task(seed=0/seed+1)`、`evaluate(bootstrap seed=0)`、`Ridge(random_state=0)` 均固定 | ✅ 通过 |
| 多重检验缺失 | 本任务为回归 rel_imp 判定，无 p 值富集 → 不适用 | ✅ N/A |
| 断言缺失 | **原代码无任何 `assert`/`isna`/`dropna`**（仅 4 处 `raise` 致命加载错误） | ⚠️ 已补（见修复） |
| 分隔符假设 | `read_csv` 4 处：2 处显式 `sep="\t"`；GSE225948 2 处原依赖逗号自动识别 | ⚠️ 已显式 `sep=","` |
| 负值/空值校验 | 无 NaN/Inf 守卫；历史上曾因 `exm1` 溢出产生 inf（已绕开 cell_ranger） | ⚠️ 补状态 score 全零守卫 |
| 重复索引 | 各 concat 前均 `var_names_make_unique()`；obs 名 make_unique | ✅ 通过 |
| 链方向/序列 | scRNA 表达分析无链方向问题 | ✅ N/A |
| LLM 幻觉包名 | 所用函数（`sc.pp.*`, `Ridge`, `MLPRegressor`）均为真实 API | ✅ 通过 |
| 版本未记录 | 原无环境锁定文件 | ⚠️ 已补 `requirements_pin.txt` |
| 中间文件未落盘 | 各步均有 print + 日志 `gate1_ts.log`；数据已落盘 | ✅ 通过 |
| 阴阳对照 | 合成 signal/null 为逻辑对照；**真实数据正对照本次补全** | ✅ 已补 |

---

## 四、生物学常识校验

- **corr=0.293**（held-out 基因表达由其他基因共表达预测）：与文献中基因共表达插值的典型 r（0.3–0.5）一致，非反常。
- **状态 score 非全零**：正对照日志显示 `DAM genes=6 / INFLAM genes=8` 进入 HVG → score 有效，排除「score 为零导致假阴性」。
- **时间轴 24h→2d→14d 质心确实有差**（方案 A 之前的质心位移任务已证实 2d≠14d）→ 数据存在真实时间生物学信号，整合无误。

---

## 五、分级审计报告

| # | 严重度 | 类型 | 位置 | 现象 / 红 flag | 对结论影响 | 状态 |
|---|--------|------|------|----------------|-----------|------|
| 1 | **阻断级** | 静默 leakage（已弃用代码仍可达） | `main.py:run_real` → `task_builder.build_real_task` | `python main.py --mode real` 仍走质心位移任务，产出**100% 虚假 PASS**（onehot(time)=答案键） | 若误用→欺诈性结论 | ✅ 已修复：`run_real` 改走 `build_heldout_task`；`build_real_task` 加 `RuntimeWarning` 弃用 |
| 2 | **重要级** | 断言缺失 | `task_builder.build_heldout_task` | 无 `state_keys` 存在性 / 状态 score 全零守卫；score 全零会静默退化为线性→假阴性 | 潜在假阴性 | ✅ 已修复：缺失列 `AssertionError`，全零 `RuntimeError` |
| 3 | **重要级** | 环境未锁定 | 全项目 | 原无依赖版本文件，pandas 3.0/numpy 2.4 较新，跨版本可能漂移 | 可复现性风险 | ✅ 已修复：`requirements_pin.txt` |
| 4 | 建议级 | 分隔符假设 | `preprocessing.py:182,192` | GSE225948 csv 依赖逗号自动识别 | 低风险（已验证可加载） | ✅ 已显式 `sep=","` |
| 5 | 建议级 | 良性警告 | `preprocessing.py:302` | `obs_names` 拼接时非唯一 `UserWarning`（`make_unique` 已在 304 行执行，无害） | 无 | ⏸ 留待（良性） |
| 6 | 建议级 | 文档过期 | `main.py` docstring / `data_acquisition.py` 注释 | 仍称 GSE174574「含 BEAM 伪时间轨迹」「真实数据走 build_real_task」 | 误导 | ✅ 已更新 `run_real` 指向 held-out |

---

## 六、🚫 投稿前必须修复清单（阻断级）

1. ~~删除/隔离 leakage 的质心位移任务入口~~ ✅ 已完成（`main.py` 已改走 `build_heldout_task`，`build_real_task` 标弃用）。
2. 论文中**不得**出现任何「质心位移任务 rel_imp=100% PASS」字样（已在 README §7 作废）。
3. Gate 1 真实结论仅以 **held-out 任务 rel_imp≈0% FAIL（真阴性）** 表述，作为「粗条件化不足 → TSC-GNN 必要」的**动机/负对照**证据。

## 七、💡 改进建议

- 将 `audit_positive_control.py` 纳入 Gate 1 标准报告流程（每次真实跑都附敏感性正对照，证明 FAIL 非假阴性）。
- `load_gse174574_raw` 拼接前显式 `obs_names_make_unique()` 消除良性警告。
- `prep_for_gate` 末尾加 `np.isfinite(adata.X).all()` 断言，拦截任何残留 NaN/Inf。
- 若将来换数据集（如人 HRA007397），复用同一套守卫 + 正对照。

---

## 八、审计结论

**正确**：任务构造无 leakage（方案 A 已通过正对照验证），评估公式对称、bootstrap 确定性、种子全固定。
**可重复**：合成与真实结果逐位复现；56/56 数据文件完整；环境已锁定；正对照证明结论非假阴性。
**可投稿状态**：Gate 1 作为「粗条件化不足」的诚实负对照 + 动机证据，**结论稳健、可写进论文**（须按第六节表述，不得引用旧 leakage PASS）。
