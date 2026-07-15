# Figure 2–7 / S1 — 数据源与可编辑格式说明

> 生成日期：2026-07-13 ｜ 作者：WorkBuddy（Repro-Guard 后续）
> 目的：记录每张图的数据来源，并把 PNG 位图转换为**可编辑的 SVG 矢量**格式。

## 一、唯一数据源脚本
所有图都由 `gate1/rewiring_study/make_figures.py` 数据驱动生成（每个 `figN()` 函数对应一张图）。
重生成命令（需 matplotlib，已装于 venv）：
```bash
cd gate1/rewiring_study
<venv>/Scripts/python.exe make_figures.py
# 同时产出 figures/*.png（稿件内嵌）和 figures/*.svg（可编辑矢量）
```

## 二、每张图的数据来源

| 图 | 函数 | 输入数据文件（同目录） | 图型 |
|---|---|---|---|
| Fig 2 Prediction benchmark | `fig2()` | **仅随机重建**（种子 20260710）；原始逐配置值未归档 | 散点 + 示意 |
| Fig 3 Temporal rewiring | `fig3()` | `rewiring_full.csv`（`q_pooled_<trans>`、`<dW_raw_><trans>`、`<dW_><trans>`） | 4 张热图 |
| Fig 4 Cross-species | `fig4()` | `human_module_enrich.json`（`per_tf[…].refs[…].OR` + 2×2 计数）、`human_bulk_gsva.json`（GSE16561 模块激活） | 森林图 + 条形图 |
| Fig 5 Drug reversal | `fig5()` | `drug_candidates_robust.csv`（`n_hits`）、`drug_perm_result.json`（置换分布） | 条形 + 直方图 + 示意 |
| Fig 6 L5 triangulation | `fig6()` | `l5_perturbation/l5_causal_direction.json`、`l5_perturbation/l5c_positive_control_result.json` | 三层卡片示意 |
| Fig 7 Evidence ladder | `fig7()` | 硬编码（来自审计 reports 文本） | 阶梯卡 |
| Fig S1 PC correction | `figS1()` | `rewiring_full.csv`（`dW_raw_<trans>` vs `dW_<trans>`） | 散点 |

## 三、可编辑格式（SVG）说明
- **格式**：`figures/figure{2,3,4,5,6,7,S1}_*.svg`，matplotlib 原生 `format="svg"` 输出。
- **可编辑性保证**：
  - `plt.rcParams["svg.fonttype"] = "none"` → 文字保留为 `<text>` 元素，可在 Inkscape / Illustrator 中直接改字。
  - 图形为 `<path>`/`<rect>` 矢量，可单独选中、改色、改形状。
  - Fig 3 热图由 `pcolormesh` 生成（**1605 个矢量单元格**），非整张栅格；仅 colorbar 渐变是 matplotlib 在 SVG 下的固定栅格行为（标准、可忽略）。
- **用途**：SVG 是真正"可编辑"的源；PNG 仍用于稿件内嵌（docx 已重建并嵌入修正后的图）。

## 四、本次同步修正（重要）
- **Fig 2 数字错误已纠正**：原 `fig2()` 硬编码写的是旧错误叙事「0 beat / 4 worse」（与稿件已修正的「4 better / 0 worse / 86 n.s.」矛盾，且方向反了）。现已改为 4 个绿点在 0 线上方、文字「4 / 90 beat · 86 n.s. · 0 worse」，与 `robustness_results.csv` 一致（REPRO_GUARD 审计，2026-07-13）。
- Fig 2A 点云为**透明重建**（原始逐配置值未归档），仅示意 4/86/0 占比，非逐点真值。

## 五、环境依赖
- venv（`C:/Users/lm962/.workbuddy/binaries/python/envs/default`）已装 `numpy 2.5.1` + `matplotlib 3.11.0`。
- ⚠️ 默认 PyPI 被公司代理 reset（10054）；须用国内镜像安装：
  `python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn matplotlib`
