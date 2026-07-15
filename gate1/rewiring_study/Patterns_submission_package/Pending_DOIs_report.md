# DOI 补齐进度报告 — TSC-GNN / 虚拟敲除 manuscript (Patterns)

生成日期：2026-07-10（续 2026-07-11 会话）
方法：Crossref REST API，分级验证（vol+page+title 高置信 / 作者+年+标题 交叉验证），
**禁止注入未经直接元数据核验的 DOI**（本次已拦截 4 个错误匹配 + 2 个会议伪 DOI）。

## 一、已填入并核验的 DOI：33 条（全部经 Crossref 直连元数据核验，title+author+year 一致）

#2, #5, #6, #7, #8, #9, #10, #11, #12, #13, #14, #20, #21, #22, #23, #24, #26, #27,
#31, #32, #33, #34, #35, #36, #41, #43, #45, #46, #51, #52, #53, #55, #56

## 二、顺带修正的 4 处引文错误（不改 DOI 指向会自相矛盾，已一并修正）

| # | 原稿错误 | 修正为（与所填 DOI 一致） |
|---|----------|---------------------------|
| 2 | Neurotherapeutics 误作 *Nat. Rev. Neurosci.* 17；标题 "repair, regeneration, remodelling" | *Neurotherapeutics* 13, 348–359；"The 3 Rs of Stroke Biology: Radial, Relayed, and Regenerative" |
| 9 | GEARS 页码 157–166 | 927–935（Nat. Biotechnol. 42） |
| 41 | scPerturb 误作 *Nat. Biotechnol.* 42, 1311–1319 | *Nature Methods* 21, 531–540 |
| 10 | scGPT 页码 1469–1480 | 1470–1480（Nat. Methods 21） |

## 三、仍待补的 23 条（按类别）

### A. 会议 / NeurIPS 论文（Crossref 无 DOI，建议改用 arXiv / OpenReview / PMLR）
- #15 Kipf & Welling 2017 GCN (ICLR) → arXiv:1609.02907
- #16 Veličković 2018 GAT (ICLR) → arXiv:1710.10903
- #17 Hamilton 2017 GraphSAGE (NeurIPS) → arXiv:1706.02216
- #18 Dwivedi 2023 "GNN: a review"（稿件标注 NeurIPS 2023，存疑）→ 经典综述为 Wu et al. 2020 *IEEE TNNLS* 10.1109/TNNLS.2020.2978386；请确认所引究竟哪篇
- #38 Zheng 2018 NO TEARS (NeurIPS) → NeurIPS 2018 / PMLR  proceedings

### B. bioRxiv 预印本（有 bioRxiv DOI，需单独查）
- #40 Bereket & Karaletsos 2022 PerturbNet → bioRxiv DOI 10.1101/2022.xxxxx（需在 bioRxiv 站查准）

### C. 期刊论文 — Crossref 标题检索无法稳定消歧（原论文未被检索排序命中，非不存在）
建议用 Zotero / EndNote 的 CrossRef 或 PubMed 直连按 vol+page 精确定位，或我改用 PubMed E-utilities 再跑一轮：
- #1 Cramer & Carrico 2008 *Nat. Rev. Neurosci.* 9, 720–731
- #3 Li 2021 *Acta Neuropathol. Commun.* 9, 152
- #4 Anrather 2024 *Nat. Immunol.* 25, 294–307
- #19 Lotfollahi 2019 scGen *Nat. Methods* 16, 1253–1261
- #25 Fancy 2011 *Ann. Neurol.* 69, 579–589
- #28 Doyle 2008 *Stroke* 39, 571–578
- #29 Jin 2010 *Prog. Neurobiol.* 90, 178–189
- #30 Barr 2010 *Stroke* 41, 2280–2285
- #37 Veličković 2023 *Nat. Rev. Phys.* 5, 343–356
- #39 Lotfollahi 2023 *Nat. Biotechnol.* 41, 1759–1770
- #42 Singh 2024 PerturbExpress *Nat. Microbiol*（稿件未给卷号，需补）
- #44 Schraivogel 2022 *Nat. Biotechnol.* 40, 1370–1378
- #47 Bhatt 2017 *Circulation* 135, 1707–1720（引文本身正确，仅 DOI 未解析出）
- #48 Li 2020 *Nat. Rev. Drug Discov.* 19, 341–359
- #49 Ahlmann-Eltze 2025 *Nat. Methods* 22, 322–331
- #50 Kartha 2023 *Nat. Genet.* 55, 1339–1350
- #54 Weider 2021 *Nat. Commun.* 12, 4240 ⚠️ 该领域高引版本为 2018 *Nat. Commun.* 9（10.1038/s41467-018-03336-3）；稿件写 2021 vol 12:4240，请核对是否同一篇/不同篇

## 四、说明与下一步
- 本次共拦截的危险错误：#16、#18 为会议论文被错配成书籍章节/综述 DOI；#47 原猜 DOI 经直连核验为 404；
  #9/#2/#41/#10 引文元数据与真实出版物不符，已修正。
- 33 条均为 vol+page+title 或 作者+年+标题 多级交叉验证通过，可放心投稿。
- 余下 23 条若需我继续，建议走 **PubMed E-utilities**（对生物医学期刊比 Crossref 标题检索更稳），
  或你用文献管理软件批量补 DOI 后我再做一轮直连核验。
