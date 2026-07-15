# -*- coding: utf-8 -*-
"""Reverse a manuscript_v7.docx back into manuscript_v7.md, faithfully preserving
the user's manual edits in the docx. Images are mapped by their order in the
document to figures/figureN_*.png (Fig1-7 inline; S1 handled by its bullet).
Display equations (centered italic LaTeX) become $$ ... $$. Inline $...$ is kept
as literal text. Superscript digits become ^N; other superscripts <sup>."""
import re, os
from docx import Document
from docx.oxml.ns import qn
from docx.text.paragraph import Paragraph
from docx.table import Table
from docx.enum.text import WD_ALIGN_PARAGRAPH

SRC = r"C:/D/workbuddy/科研/虚拟敲除/gate1/rewiring_study/manuscript_v7.docx"
OUT = r"C:/D/workbuddy/科研/虚拟敲除/gate1/rewiring_study/manuscript_v7.md"

FIGFILES = {
    1: "figure1_tsc_gnn_conceptual_framework_v4.png",
    2: "figure2_prediction_benchmark.png",
    3: "figure3_rewiring_heatmaps.png",
    4: "figure4_cross_species.png",
    5: "figure5_drug_reversal.png",
    6: "figure6_L5_triangulation.png",
    7: "figure7_evidence_ladder.png",
    8: "cellchat_rewiring_heatmap.png",
}

doc = Document(SRC)
body = doc.element.body

# ---- pre-scan blocks ----
blocks = []
for child in body.iterchildren():
    if child.tag == qn('w:p'):
        p = Paragraph(child, doc)
        drawing = child.find('.//' + qn('w:drawing')) is not None
        blocks.append(('p', p, drawing))
    elif child.tag == qn('w:tbl'):
        blocks.append(('tbl', Table(child, doc), False))


def inline(p):
    segs = []
    for r in p.runs:
        t = r.text
        if not t:
            continue
        bold = bool(r.bold)
        ital = bool(r.italic)
        if r.font.superscript:
            t = ("^" + t) if re.fullmatch(r"\d+", t) else f"<sup>{t}</sup>"
        elif r.font.subscript:
            t = f"<sub>{t}</sub>"
        segs.append((bold, ital, t))
    merged = []
    for b, it, t in segs:
        if merged and merged[-1][0] == b and merged[-1][1] == it:
            merged[-1][2] += t
        else:
            merged.append([b, it, t])
    res = []
    for b, it, t in merged:
        if b and it:
            t = f"***{t}***"
        elif b:
            t = f"**{t}**"
        elif it:
            t = f"*{t}*"
        res.append(t)
    return "".join(res)


def render_table(table):
    rows = table.rows
    if not rows:
        return ""
    data = []
    for row in rows:
        cells = []
        for cell in row.cells:
            txt = " ".join(pp.text for pp in cell.paragraphs if pp.text.strip())
            cells.append(txt)
        data.append(cells)
    lines = []
    ncol = len(data[0])
    lines.append("| " + " | ".join(data[0]) + " |")
    lines.append("| " + " | ".join(["---"] * ncol) + " |")
    for r in data[1:]:
        lines.append("| " + " | ".join(r) + " |")
    return "\n".join(lines)


out = []
title_done = False
in_supp = False
skip_supp = False


def emit(line):
    if line == "":
        if out and out[-1] != "":
            out.append("")
    else:
        out.append(line)


for i, blk in enumerate(blocks):
    kind = blk[0]
    if kind == 'tbl':
        emit(render_table(blk[1]))
        emit("")
        continue
    p = blk[1]
    drawing = blk[2]
    text = p.text
    style = (p.style.name if p.style else "") or ""
    is_center = (p.alignment == WD_ALIGN_PARAGRAPH.CENTER)

    if drawing:
        nxt = ""
        for j in range(i + 1, min(i + 6, len(blocks))):
            if blocks[j][0] == 'p' and blocks[j][1].text.strip():
                nxt = blocks[j][1].text.strip()
                break
        if re.search(r"Supplementary Figure\s+S\d+", nxt):
            continue  # S1 embedded via its bullet; skip inline
        mfig = re.match(r"Figure\s+(\d+)", nxt)
        if mfig:
            N = int(mfig.group(1))
            if N in FIGFILES:
                emit(f"![Figure {N}](figures/{FIGFILES[N]})")
                emit("")
                continue
        # fallback ordinal
        img_ord = sum(1 for b in blocks[:i] if b[0] == 'p' and b[2]) + 1
        if img_ord in FIGFILES:
            emit(f"![Figure {img_ord}](figures/{FIGFILES[img_ord]})")
            emit("")
        continue

    if not text.strip():
        emit("")
        continue

    # title
    sz = None
    for r in p.runs:
        if r.font.size is not None:
            sz = r.font.size.pt
            break
    if (not title_done) and style == "Normal" and is_center and sz == 18 and any(r.bold for r in p.runs):
        emit("# " + text)
        title_done = True
        continue

    if style.startswith("Heading"):
        lvl = int(re.sub(r"\D", "", style) or 1)
        t = text.strip()
        if t == "Supplementary Materials":
            in_supp = True
        if t == "References":
            in_supp = False
        if lvl == 1:
            emit(("# " if in_supp else "## ") + t)
        elif lvl == 2:
            emit(("## " if in_supp else "### ") + t)
        elif lvl == 3:
            emit("#### " + t)
        continue

    if style == "List Bullet":
        emit("- " + inline(p))
        continue
    if style == "List Number":
        emit("1. " + inline(p))
        continue

    # display equation: centered + italic + math chars, no inline $
    if is_center and any(r.italic for r in p.runs) and "$" not in text and re.search(r"[\\{}_^=]", text):
        emit("$$ " + text.strip() + " $$")
        emit("")
        continue

    emit(inline(p))

# collapse trailing blank lines
while out and out[-1] == "":
    out.pop()

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
print("wrote", OUT, "| lines:", len(out))
