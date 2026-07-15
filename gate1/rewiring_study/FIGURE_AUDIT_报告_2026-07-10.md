# Figure-生成流程审计报告（research-bug-audit · 模式 A）

- **审计对象**：`gate1/rewiring_study/make_figures.py`（生成 Fig 2–7 + S1）
- **审计日期**：2026-07-10
- **模式**：模式 A（既往流程审计，找已有 bug + 复现性核对）
- **运行环境**：conda `bbb_gnn`（matplotlib 3.10.9 / numpy 2.2.6）；managed venv 无 matplotlib

---

## 一、结论速览

**复现性：通过。** 在相同环境下重跑两次，7 张 PNG 与原始生成版（20:48）**逐字节一致（md5 相同）**，运行无报错/警告。Fig 2 用 `np.random.default_rng(20260710)` 固定种子，确定性良好。

**数值正确性：核心量化面板正确。** Fig 4 的 OR（27.0 / 16.5 / 4.6）与 `human_module_enrich.json` 中预计算值逐一吻合；Fig 6 的 L5 数值（Sox10 46/412·OR↓1.81、Cebpb 149/404·OR↓1.20、MYC 19/332、BCL11A 29/332、GATA1 187/332）与 `l5_*.json` 逐一吻合；Fig 5B 置换（observed=4, emp-p=0.33, n=20）与 `drug_perm_result.json` 吻合；Fig 3 显著边数（8/28/19/36）与正文 Table 2 吻合。

**需修/需核对：4 项（见下表 + 文末清单）。** 无会导致错误结论的"静默数值 bug"，但有 1 项图形缺陷（Fig 5A 平条）和 2 项图文一致性问题，以及多处复现性加固点。

---

## 二、审计发现表

| 序号 | 严重度 | 类型 | 位置（文件:行） | 现象 / 红 flag | 对结论的潜在影响 | 修复建议 |
|------|--------|------|----------------|--------------|------------------|----------|
| 1 | 重要 | 图形缺陷 / 可能误导 | `make_figures.py:255-273` | Fig 5A 取 CSV **前 15 行**（按文件 `rank`）画 `best_score`；但 CSV 的 `rank` 并非按 `best_score` 排序——rank 1–16 全部 =0.125、17–26=0.15、27–30=0.175。故前 15 条**全部平齐在 0.125**，条形图退化成 15 根等高柱；有区分度的分数（0.15/0.175）在 rank 17–30 被排除。 | "Top robust L1000 compounds" 看似有排名，实际展示子集无差异，审稿人易视为"图坏了"或暗示不存在的梯度。4 个文献支持药（红柱）信息仍在，但图整体信息量低。 | 按 `best_score` 降序取 top15（或改画有方差的 `n_hits` 1–6 / `mean_score` 0.125–0.14）；或图注明确"分数为分层阈值，前 15 名并列"。 |
| 2 | 重要 | 复现性 / 硬编码未读源 | `make_figures.py:190-203`（Fig4A）、`316-329`（Fig6）、`367-378`（Fig7） | OR 面板硬编码 `(k,n,nt)` 与 `N=18564`；L5 数值全部硬编码。声明的源文件 `human_module_enrich.json` / `l5_*.json` **存在但未被读取**（L5 json 实际在 `l5_perturbation/` 子目录）。 | **当前数值已逐一核验正确**，但图形不随分析重跑自动更新——若上游富集/L5 分析重算，图会静默偏离源数据。属"将来复现"风险。 | Fig4A 改为从 `human_module_enrich.json["per_tf"][TF]["refs"][ref]` 读 `k/n/OR`，从 `["N_universe"]` 读 N；Fig6 从 `l5_causal_direction.json` / `l5c_positive_control_result.json` 读数值；并加一条交叉校验断言。 |
| 3 | 重要 | 图文一致性 | 正文 L206 vs Fig S1 数据 | 正文称"77–83 % of edges preserved their direction after PC composition correction"；但 Fig S1 所用边集（pooled q<0.05）实测方向保持率：**汇总 73.6%**，分过渡 75.0% / 75.0% / 78.9% / **69.4%**（sham→14d 明显低于 77%）。 | 审稿人对照 Fig S1 与正文会发现数字对不上；疑似正文语句来自更早版本的分析边集。 | 核对 77–83% 的精确算法/边集；若与图不一致，把正文改为实际的分过渡率（或注明所用边集定义），使图文自洽。 |
| 4 | 重要 | 诚实报告（视觉-文字自相矛盾） | `make_figures.py:76-98` | Fig 2A 用种子随机散点：86 个"无差异"点取 N(0, 1.4)，约一半落在 0 线之上；但面板文字断言"0/90 beat linear · 86/90 not significant"。视觉印象（许多点在 0 线上方）与文字声明（0/90 胜出）冲突。 | 稿件以"诚实边界"为核心卖点，审稿人一眼可见图内自相矛盾，削弱可信度。已标注"transparent reconstruction"，但视觉仍误导。 | 把"无差异"云中心下移（如 N(−0.5, 1.0)）使少数越线；或在图中加注"纵轴位置为示意抖动，唯 0/86/4 计数为事实"。 |
| 5 | 重要 | 复现性 / 环境未锁定 | 全局 | 脚本仅在 `bbb_gnn`（matplotlib 3.10.9 / numpy 2.2.6）可跑；managed venv 无 matplotlib；脚本未记录版本、未写日志。 | 在别处运行会产生不同字节的 PNG（字体度量差异）；复现脆弱。 | 脚本顶部注释固定环境（bbb_gnn + 版本）；可选在 `__main__` 打印 `matplotlib.__version__`/`numpy.__version__` 到日志。 |
| 6 | 建议 | 鲁棒性（潜在崩溃） | `make_figures.py:177-178`（Fig3）、`426`（S1） | Fig3 `fig.colorbar(im)` 的 `im` 仅在 `if el:` 内定义——若某次运行**所有过渡**均无显著边，`im` 未定义 → NameError。Fig S1 `min(raw+corr)` 在空列表上 → ValueError。 | 当前数据（8/28/19/36 边）不触发；但空输入会崩。 | 对 `im` 与 `raw/corr` 加空值保护（无显著边时跳过色条/标注"无显著边"）。 |
| 7 | 建议 | 静默失败预防 | 全局 | 无断言校验输入文件/列存在、无 `gs["modules"]` 含预期模块的校验。若 CSV 分隔符/列名变更 → 热图变空、模块被静默丢弃。 | 分隔符风险已排除（实测为逗号），但无防护，将来换数据易静默出错。 | 开头加断言：各输入文件存在、关键列存在、`human_bulk_gsva.json` 含 SOX10/CEBPB/GATA2/PAX6。 |
| 8 | 建议 | 数据卫生 | `drug_candidates_robust.csv` 表头 | `rank` 列带 BOM：`\ufeffrank`。对脚本实际使用的 `drug`/`best_score` 无影响，但一旦改用 `row["rank"]` 会 KeyError。 | 潜伏地雷。 | 读取用 `utf-8-sig`，或修正源 CSV 去 BOM。 |
| 9 | 建议 | 代码整洁 | `make_figures.py:15` | `import random` 未使用（图用 `np.random.default_rng`）。 | 无害。 | 删除未用导入。 |
| 10 | 建议 | 鲁棒性 | `make_figures.py:241` | Fig4B `axB.set_ylim(0, 0.62)` 硬编码；当前最大 stroke_mean=0.529 安全，但数据增长会裁切。 | 当前安全。 | 从数据计算 ylim（留 10% 余量）。 |

---

## 三、🚫 投稿前必须修复 / 核对清单（阻断级 + 重要级）

1. **Fig 5A 平条缺陷（#1）**：改为按分数排序或改画 `n_hits`，使图有信息量、不自相矛盾。
2. **Fig S1 ↔ 正文 L206 数字不一致（#3）**：核对 77–83% 算法，使图文自洽（建议正文改成分过渡实测率 75/75/79/69% 或注明边集）。
3. **Fig 4A / Fig 6 硬编码改读源（#2）**：为"将来复现"计，改为从 JSON 读取并在脚本内做交叉校验断言（当前数值已正确，不改不影响本期结论，但强烈建议做）。
4. **Fig 2A 视觉-文字冲突（#4）**：下移"无差异"云或加注，消除图内自相矛盾。

> 说明：严格意义上**未发现"静默数值错误"型阻断 bug**（FDR 已在 Fig3/S1 用 pooled q<0.05 正确施加；Fig2 种子已固定；各量化面板数值已与源 JSON 逐一核对正确）。上述 4 项属于**图形质量 / 图文一致性 / 复现加固**，但因稿件以"诚实边界"为卖点，建议在投稿前处理。

---

## 四、💡 改进建议（非阻断）

- #5 固定运行环境（bbb_gnn + 版本注释 / 版本日志）。
- #6 空显著边集时的崩溃保护。
- #7 输入文件/列/模块的断言校验。
- #8 去 BOM（读取用 utf-8-sig）。
- #9 删 `import random`。
- #10 Fig4B ylim 由数据计算。

---

## 五、复现性验证记录（已通过）

| 检查项 | 结果 |
|--------|------|
| 输入文件齐全（rewiring_full.csv / human_bulk_gsva.json / drug_candidates_robust.csv / drug_perm_result.json / human_module_enrich.json / l5_perturbation/l5_*.json） | ✅ 全部存在 |
| CSV 分隔符 | ✅ 逗号（非 tab），关键列 `dW_/dW_raw_/q_pooled_` 齐全 |
| 重跑无报错/警告 | ✅ 7 张 `saved` + `ALL FIGURES DONE` |
| 确定性（RUN1 vs RUN2） | ✅ **逐字节一致** |
| 与原始 20:48 PNG 一致 | ✅ **逐字节一致** |
| Fig 4 OR（27.0/16.5/4.6） | ✅ 与 JSON 预计算值吻合 |
| Fig 6 L5 数值 | ✅ 与 `l5_*.json` 吻合 |
| Fig 5B 置换（4 / 0.33 / n=20） | ✅ 与 JSON 吻合 |
| Fig 3 显著边数（8/28/19/36） | ✅ 与 Table 2 吻合 |
| Fig S1 r(raw, corrected)=0.02 | ✅ 与正文"方向保持 77–83%、部分翻转"框架自洽（原始与校正近乎正交，说明校正移除主导混杂，非图 bug） |

---

## 六、生物学常识校验

- Fig 3：ΔW 热图来自真实 rewiring 输出，边数匹配正文；合理。
- Fig 4：SOX10→髓鞘（OR 27）、CEBPB/GATA2→神经炎症（OR 16.5/4.6）方向为正、与已知生物学一致；人血激活 q≈1.1e-8 且 PAX6 正确标为负对照——合理。
- Fig 5：命中含 HDAC 抑制剂（vorinostat/trichostatin A）与他汀（mevastatin/rosuvastatin），有文献支撑；置换 p=0.33（n.s.）如实展示——合理。
- Fig 6/7：三层三角验证与 L1–L5 阶梯与报告一致——合理。
- Fig S1：原始 vs 校正 ΔW 近乎正交（r=0.02），与正文"校正移除组成混杂、更正估计更可靠"的论述一致，非错误。

---

## 七、建议的后续动作（待作者确认）

A. 应用 #1/#4 图形修复并重生成 7 图（会改变 PNG 字节，但不改结论）。
B. 应用 #2 把 Fig4A/Fig6 改为读 JSON + 交叉校验（提升复现性，数值不变）。
C. 核对 #3 正文 L206 数字并修订（图文自洽）。
D. 应用 #5–#10 加固（环境锁定 + 断言 + 整洁）。

以上 A–D 是否执行、按哪种顺序，请你定夺；B/C/D 不改本期科学结论，A 仅改可视化呈现。
