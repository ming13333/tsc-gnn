# -*- coding: utf-8 -*-
"""Generate L5 three-layer triangulation as a ProcessOn .pos file."""
import sys, os
sys.path.insert(0, r"C:\Users\MI\.workbuddy\skills\generate-pos-flowchart\scripts")
from generate_pos import PosBuilder

B = PosBuilder()

# Title + legend
B.add_text("L5 三层三角验证：因果支持仅成立于恰当生物学语境", 80, 18, w=640, h=44, font_size=24)
B.add_text("绿色 = 阳性证据    琥珀 = 语境边界 null", 80, 66, w=640, h=24, font_size=14)

# Three evidence-layer cards (top row)
L5a = B.add_process(
    "L5a 原生谱系 bulk KO\n\nSox10 cKO / Cebpb het-KO\n靶程序在自身 KO 下显著下调\nSox10 rank 46/412 (top 11%)\nCebpb rank 149/404 (top 37%)\n\n+ 阳性 (native-lineage)",
    40, 110, w=240, h=190, fill_color="225,245,225", font_color="27,94,32")

L5b = B.add_process(
    "L5b LINCS 过表达 OE\n\nGATA2 OE vs DoRothEA 靶程序\nmimicker, 方向一致\nrank 3 / 33,782 (top 0.01%)\np = 1.4e-05, 强自特异性\n\n+ 阳性 (OE, 正确方向)",
    280, 110, w=240, h=190, fill_color="225,245,225", font_color="27,94,32")

L5c = B.add_process(
    "L5c K562 sc-CRISPRi\n\nReplogle 2022 genome-scale\nCEBPB self-Z = -0.43 → 靶集不动 (p=0.84)\nGATA2 self-Z = -0.19 → 不动 (p=0.18)\nSOX10 位点不在 K562 基因空间\n\no 语境边界 null (off-context)",
    520, 110, w=240, h=190, fill_color="255,248,225", font_color="180,120,10")

# Center conclusion card
CTR = B.add_process(
    "三模态独立收敛 →\n\n因果支持仅成立于\n恰当生物学语境\n\nL5 = causal support\n≠ validation",
    250, 350, w=300, h=160, fill_color="150,60,200", font_color="255,255,255")

# Positive-control footnote strip
FN = B.add_process(
    "K562 内部正对照 (方法有效功率校准)：\nMYC rank 19/332 (p=3.1e-3)、BCL11A rank 29/332 (p=7.5e-3) → 靶集显著下移\n⇒ L5c 的 null 是真语境边界，非流程失效；\n但 GATA1 (on-target 敲低成功) 靶集仍不动 → 通用 regulon 仅粗代理",
    40, 560, w=720, h=110, fill_color="255,250,235", font_color="120,80,0")

# Converging links: three layers -> center
B.add_linker(L5a, CTR, from_anchor=1, to_anchor=0)
B.add_linker(L5b, CTR, from_anchor=1, to_anchor=0)
B.add_linker(L5c, CTR, from_anchor=1, to_anchor=0)
# center -> positive-control (supporting relationship)
B.add_linker(CTR, FN, from_anchor=1, to_anchor=0)

data = B.build("L5 三层三角验证 (A+C 模式)")
out = r"C:\D 盘\科研\虚拟敲除\gate1\rewiring_study\l5_perturbation\L5_triangulation.pos"
B.save(data, out)
