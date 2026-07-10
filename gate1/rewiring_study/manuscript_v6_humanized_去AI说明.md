# 去 AI 痕迹说明（manuscript_v6 → manuscript_v6_humanized）

应用 `academic-de-ai-workflow`（humanizer + humanizer-nature-style + nature-polishing）对 Patterns 格式的 v6 稿件做了多轮去 AI 处理。**原稿未覆盖**，输出为新文件 `manuscript_v6_humanized.md`。所有科学内容（数据、基因名、p 值、OR、ρ、Jaccard、rank、参考文献编号、图表、公式）均保持不变。

## 处理范围与保留项
- **受保护、未改动**：`## References` 整段（Cell Press 式，DOI 仍 `[DOI: pending]`）；所有表格行（含单元格内的 `—` 占位符）；`![Figure ...]` 图片标签；`[Figure X about here — ...]` / `[Table X about here — ...]` / `[Author note — ...]` 等占位符（保留其 em dash）；`$$...$$` 与 `$...$` 数学公式；Terminology Ledger、Highlights、Keywords。
- **处理对象**：正文散文（Introduction / Methods / Results / Discussion / Conclusion / 图注 / Bigger picture / eTOC / blockquote 披露框）。

## Round 1 词汇/句法层（humanizer）
- em dash（—）：正文 106 处 → 1 处（仅 Table 10 单元格 `rank 304/332 — locus absent` 作为表格占位符保留）。其余改为逗号/句号。
- 删除填充词/AI 高频副词：`Moreover`/`Furthermore`/`Additionally`/`Notably`/`Importantly`/`additionally` 全部移除。
- `we argue` ×2：移除（含修复遗留逗号，"that, benefits" → "that benefits"）。
- `serves as` → `is`；`highlights`(动词) → `shows`。
- `closed the loop` → `links to`（×2）；`a closed loop` → `a complete path`。
- `most generalisable` → `central`；`where this work begins, not where it fails` → `the starting point of this work`。
- 防御性框架：`deliberate design choice, not a limitation` → `intentional`；`is deliberate and is itself a key finding` → `is itself a key finding`；`openly`/`transparently` 移除（"report openly" → "report"）。

## Round 2 结构/逻辑层（humanizer-nature-style）
- 移除元叙事与"工具自述"：`recovery engine` / `prediction engine` 标签重写为"recover, not predict" / "recovery-oriented method rather than a prediction method"。
- 修复 em dash→逗号造成的逗号粘连（Methods 2.2 一句拆为两句）。
- 改写 3.1 开头元叙事（"rest on a scoping result that first fixes..." → "build on a scoping result that establishes..."）与 3.2 的 Q&A 框架（"This raises an immediate question... the framework is the instrument..." → "The graph earns its place not by improving prediction but by enabling edge-level recovery..."）。
- `interpretability-first`：4.1 处保留语义后重写为 "a virtual-perturbation framework built around interpretability"；3.8.2 处改为 "an interpretable-structure framework"。

## Round 3–4 润色 + 术语一致
- `progressively strengthen confidence` → `build confidence`（图 1 图注）。
- 术语一致：DoRothEA / TSC-GNN / GRN 等全程统一；引用句动词多样化；Discussion 各段均绑定自有数据（L1–L5）。

## 验证（grep 核对）
- **已消除**：recovery engine、prediction engine、we argue、deliberate design choice、most generalisable、closed the loop、where this work starts、Moreover/Furthermore/Additionally/Notably/Importantly、正文散落 em dash。
- **完整保留**：p = 0.022 / 0.0005、OR = 27、ρ = 0.48–0.55、Jaccard ≈ 0、ΔW = +0.51、rank 3/33,782、rank 46/412、SOX10/CEBPB/GATA2、<sup> 引用、Table S1–S8、Fig. 7、0/90、n_perm。

## 需作者确认的判断点
1. **Discussion "Anticipated reviewer questions"（Q1–Q3）保留**：这是方法学论文常见的主动回应审稿人结构，非生成式 AI 的会话填充。如需进一步去除，可改为连续段落。
2. **Table 10 单元格 `—` 保留**：属表格占位符，非散文 em dash。
3. **参考文献 DOI 仍为 `[DOI: pending]`**：属格式遗留（非 AI 痕迹），投稿前须补真实 DOI 与期刊标准缩写。
4. 全文主文约 11.5k 词，Patterns 无硬上限但典型 5–7k；如需更贴合可再压（本次未动长度以保证据深度）。
