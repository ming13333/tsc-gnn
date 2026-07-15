# 可复现性与可靠性审计（AI Research Reproducibility Guard, Reviewer Mode）

**日期**：2026-07-13 ｜ **对象**：`manuscript_v6_humanized.md` + `gate1/` 全部分析脚本与数据
**方法**：独立复核原始脚本与 CSV，不采信任何摘要/记忆中的结论。Reviewer-mode，不夸赞、只列问题。

---

## 0. 总判定

| 维度 | 状态 | 说明 |
|---|---|---|
| **可靠性（结论是否站得住）** | 🟢 **致命项已修复** | 核心 scoping 数字 0→4/90 已精核并改正；主轴（图无稳定预测优势→主打可解释性）成立 |
| **可复现性（结果能否复现）** | 🟡 **部分达标** | 种子已固定、图可逐字节复现；但环境锁文件互相矛盾、数据来源未记版本、manifest 声明过宽、LINCS 未定版（B1–B6 待办）|

**一句话**：数据本身大体扎实且多有诚实披露（L4/L5c 尤其）。**核心 scoping 结论原在稿中写成"0 better"是数字错误（正确为 4/90，全来自单一种子的组合任务），已修正**；本 agent 初次审计误把整表"24 better"当反驳，亦属口径混淆，一并更正。可靠性风险已消除，剩余为 B1–B6 可复现性工程项。

---

## A. 可靠性发现

### A1. 🔴【致命→已修复】中心结论数字写错：应为「4/90」而非「0」

> **2026-07-13 二次精核后更正（含更正本 agent 自己在 A1 初稿中的过度表述）。**

**逐行重核（按报告自定义的口径分类，而非笼统数整表）**：

| 口径 | 定义 | better / worse / n.s. | 判定 |
|---|---|---|---|
| **报告"90 点"主基准** | 图变体 × **线性 ridge 读头**（single 5图×10=50 + combo 4图×10=40 = 90）| **4 / 0 / 86** | ✅ 报告 §五正确 |
| 稿件 §3.1/§4.1/Summary | 同上 90 config | 写成 **0 / – / 86** | ❌ 数字错，应为 4 |
| 整表 192 行（含消融） | +`+deg`(ridge no-P) 96 行 + kernel 探针 6 行 | 24 / 17 / 151 | ⚠️ 不能用来反驳"90 config" |

- **正确数字 = 4**（不是稿件的 0，也不是本报告初稿误引的 24）。那 4 个 better **全部来自 seed1 的 combo 任务**（gnn_knn/gnn_doro/gnn_rand/0-hop，rel≈+18~19%，CI 不含 0）；其余 9 个种子的 combo 全 n.s.，single 50/50 全 n.s.，**worse=0**。
- **本报告初稿的错误**：把整表的 "24 better" 当作对稿件"90-config = 0 better"的反驳。但那 24 里 **20 个来自 `+deg`（抽掉扰动 ID 通道）的旁支消融**，不属于稿件定义的 90-config 主基准。用它反驳，与稿件写"0"是**同一类口径混淆**——特此更正。
- `+deg` 消融本身确有信号（combo 上 20/40 `better`：抽掉扰动 ID 后图传播确实补偿了预测），属**有价值的附带发现**，但它是消融而非主基准，不应进入"90 config"标题，也不足以推翻主轴。

**对主轴的影响（修正后）**：主轴**不崩**。把 0 改为"4/90（且全来自单一种子的组合任务）、其余 86 无差异、无一劣于线性"后，"图无**稳定**预测优势 → 主打边级可解释性"依然成立，且比写 0 更诚实、更抗审稿。

**已执行（有界事实修正，非改叙事）**：稿件 6 处已由 0/"no gain" 改为 "no consistent gain / 4-of-90"（Highlights L19、eTOC 邻近、Summary L30、§2.10 L148、§3.1 L156、§4.1 L338、Bigger picture L38）。docx 已重建、投稿包已同步。

### A2. 🟠【严重】错误结论已被工作记忆固化（AI 静默错误）

工作记忆 `MEMORY.md` 与本轮对话摘要均记载：
> "图优于线性 0 例、无差异 86 例 →『固定图+线性读头图未带来预测优势』为强结论"

这正是 A1 的错误被当作事实在多个会话间传递。**AI 不会自我怀疑**，把 CSV 上方的错误解读一路带成了「强结论」。印证了 repro-guard 的核心警告：agent 输出必须人工独立复核，绝不可让 agent 自审。

### A3. 🟡【轻微】L4「robust hits」措辞略过声称

- L4 表（L185）写 "Robust hits vorinostat / trichostatin A / Mevastatin / Rosuvastatin"，但同表已标注 permutation **p = 0.33 (n.s.)**。
- 正文（L250–254）其实**很诚实**：明确 "empirical p ≈ 0.3 … do not claim significant enrichment"，并有独立「Author note」限制框（hypothesis-generation, not efficacy prediction）。
- **建议**：把表中 "Robust hits" 改为 "Candidate / top-50 hits (p = 0.33, n.s.)"，与正文一致即可。不构成致命问题。

### A4. 🟢【优点】L5c 阴性结果处理诚实

- Replogle K562 CRISPRi 中三 TF 靶程序在自身敲低下**无显著下移**，被框架为「语境边界 null」而非管线失败。
- 有专用正对照脚本 `l5c_positive_control.py` 与 docstring "NOT validation"。
- 这是合规的、可接受的，保留。

### A5. 🟡【注意】A+C「把表面阴性转正面证据」属事后拟合

- 记忆记录 "梯度即结果=框架预测，把表面阴性转正面证据"。
- 作为**假设**可接受，但不得表述为「验证」。当前稿件 Scope disclosure（L260）已声明 "directional causal support, explicitly not edge-level validation"，基本合规；仅需在讨论中避免把「语境门控恰为框架预测」写成因果证据。

---

## B. 可复现性发现

### B1. 🟠【严重】环境锁文件互相矛盾

| 文件 | python | numpy | pandas | 备注 |
|---|---|---|---|---|
| `tsc-gnn-repo/environment.yml` | 3.11 | 2.2.6 | 2.3.3 | scikit-learn 1.5.*，torch **未定版**，decoupler/mygene/requests 未定版 |
| `gate1/requirements_pin.txt` | — | **2.4.6** | **3.0.3** | scikit-learn 1.9.0、scanpy 1.12.2——版本疑似超前/不可装 |
| `run_manifest.txt`（实测） | **3.10.20** | 2.2.6 | 2.3.3 | 与 env.yml 的 3.11 不一致 |

- `requirements_pin.txt` 的 numpy 2.4.6 / pandas 3.0.3 在当前时点**很可能不存在或不可装**，该锁文件不可信。
- torch（GEARS 比对必需）在 environment.yml 中**未定版**。
- **必须做**：统一为一个可安装的 `environment.yml`（含精确 torch / decoupler / scanpy / mygene 版本），删除或修正 `requirements_pin.txt`。

### B2. 🟠【严重】数据来源版本/日期未记录（provenance 缺失）

- GEO 下载脚本（`data_acquisition.py` 等）**未记录 accession 版本或下载日期**。重拉 GSE174574/225948 可能拿到不同元数据/注释。
- DoRothEA 仅本地 `.tsv`，**无版本/释放号**（且存在未使用的 `mouse_dorothea_regulon_v1.0.tsv` 99B 占位）。
- Replogle h5ad、GSE269122/273163 均为本地文件，**无版本字符串**。
- **必须做**：在每个数据集旁写 `DATA_PROVENANCE.txt`（accession + version/date + checksum），或在 manifest 中登记。

### B3. 🟠【严重】LINCS 查询未定版 + 单次实时 API

- `drug_reversal_permutation.py` L50：`"db-version": "latest"` → **未冻结**，重跑可能换库。
- L252 自述 "single-shot, live API" → 该置换**无法确定性复现**。
- N_PERM=20 偏低（最小可达 p≈0.048）。
- **建议**：固定 `db-version` 具体号；把 L1000CDS2 结果落盘为 JSON 归档，置换从落盘文件重算，避免依赖实时 API。

### B4. 🟡【中等】「All runs emit a manifest」声明过宽

- 仅 `run_rewiring.py` 落 `run_manifest.txt`（含命令/版本/seed/缓存 sha256，质量尚可）。
- `robustness_study.py`、`drug_reversal_permutation.py`、L5 系列脚本**均无 manifest 输出**。
- **必须做**：稿件 §2.9（L144）"All runs emit a manifest" 应改为 "The rewiring analysis emits a manifest …"；或给其余脚本补齐 manifest。

### B5. 🟡【中等】run_manifest 缓存路径已失效

- `run_manifest.txt` L20：`cache_path: C:/D 盘/科研/虚拟敲除/.../rewiring_doro_cache.npz`——**旧机器路径，现机不存在**（`C:/D/workbuddy/...`）。
- 记录的 `cache_sha256` 指向一个**当前不可达**的文件 → 该缓存的复现性声明实际落空。
- **必须做**：更新路径并把 `.npz` 一并归档到仓库/补充材料，或重建后重算 hash。

### B6. 🟡【中等】Phase-2 GNN 核心未固定种子

- `gate1/tsc_gnn/model.py`、`train_eval.py`、`run_phase2.py` 中**未发现 seed 设定** → 报称的 73.8% Phase-2 结果**不受种子控制**。
- **必须做**：补 `torch.manual_seed` 并报告多种子区间，或明确标注其为示意性/非种子锁定结果。

### B7. 🟢【优点】图件可复现

- `FIGURE_AUDIT_报告_2026-07-10.md`：同环境重跑 2 次，7 张 PNG **md5 逐字节一致**。确定性良好，保留。

### B8. 🟢【优点】存在审计痕迹

- `gate1/README.md`（明令禁止引用旧 100% PASS 泄漏结果）、`AUDIT_REPORT_2026-07-09.md`、多个 run log、`tsc-gnn-repo` git 历史（127 文件、20+ 提交）。
- 部分满足"Archive artifact"要求，但缺独立的「prompts→code→run→human-verify」连贯 log。

---

## C. 投稿前可复现性清单（checklist §5）

| Item | 状态 | 行动 |
|---|---|---|
| 全部引文 DOI 已核验 | 🟡 | 16 条已修；**7 条待你决策**（含 4 条疑似虚构，见 DOI 核验报告） |
| 全部数据 raw source 已存 | 🔴 | 缺 GEO/DoRothEA/LINCS 版本与路径（B2/B5） |
| 代码 git 版本化/冻结 | 🟢 | git 历史齐全 |
| 环境锁版本固定 | 🔴 | 锁文件矛盾（B1） |
| 随机种子固定并记录 | 🟡 | 部分固定；Phase-2 缺失（B6） |
| 统计方法专家复核 | 🟡 | 核心结论数字已精核修正（A1，0→4/90）；仍建议统计人再审 CI/置换 |
| 模型参数记录 | 🟡 | 部分 |
| 图由代码生成 | 🟢 | md5 一致 |
| 关键结果独立验证（第二法） | 🟢 | A1 已由脚本按口径重算复核（4/0/86） |
| AI 生成文本人工批准 | 🟡 | 需你最终逐句审 |

---

## D. 人工必须执行的核查（DRIL：Human verifies）

1. **【已修复，请复核】A1 数字更正**：正确数字 = **4/90**（base-graph×ridge 主基准，全来自 seed1 combo；worse=0）。整表的 24 better 含 +deg 消融不可用作反驳。稿件 6 处 0→"4/90/no consistent" 已改、docx 已重建。**请你核对 `robustness_results.csv` 的 `sig` 列确认 4 无误；若你更希望保留"仅在最保守设定下无优势"的更强表述或反向重定位主轴，再告诉我。**
2. **A3 措辞**：确认把 L4 表 "Robust hits" → "Candidate hits (p=0.33, n.s.)"。
3. **B1 环境**：统一 `environment.yml`，删/修 `requirements_pin.txt`，补 torch 定版。
4. **B2/B5 provenance**：补 `DATA_PROVENANCE.txt` + 更新缓存路径与 hash。
5. **B3 LINCS**：固定 `db-version`，落盘 API 结果。
6. **B4/B6 manifest/seed**：收窄 §2.9 声明或补齐 manifest；Phase-2 补种子。
7. **7 条存疑引文**：审阅 `DOI核验与引文修正报告_2026-07-13.md` 后逐条拍板（含 4 条疑似虚构是否删除）。

> 原则（repro-guard 核心）：**研究结论的真伪归您，agent 只做有界执行**。本审计中 A1 类矛盾若由 agent 自行"修正叙事"会再次引入静默错误——故所有结论性改动须您拍板。

---

## E. 本 agent 此前犯错的透明说明

- **错误一**（更早会话）：将 "0 better / 86 n.s." 当作已验证事实写入 `MEMORY.md` 与稿件，**未回读 CSV 的 `sig` 列**。
- **错误二**（本审计初稿）：纠正错误一时用力过猛，把**整表 24 better**（含 20 个 `+deg` 消融行）当作对"90-config"主基准的反驳，声称"图在组合扰动上确实更好、反而支持图方法"——这是**第二次口径混淆**，同样未按报告自定义的子集分类。
- **最终精核**：报告"90 点"主基准 = base-graph×ridge = **4 better / 0 worse / 86 n.s.**；`+deg`/kernel 属消融，另计。稿件应写 4。两次错误均已更正。
- **教训**：任何"强/致命结论"落盘前，必须(1)回读原始数据，(2)**严格对齐口径定义**（哪些行进入分母），(3)避免为"翻案"而反向夸大。
