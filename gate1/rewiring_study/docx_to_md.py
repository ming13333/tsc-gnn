# -*- coding: utf-8 -*-
"""Extract a faithful markdown mirror from manuscript_v7.docx.

The docx is now the canonical (manually edited) artifact, so we regenerate the
.md so the repo carries a diffable text version. We preserve:
  - Title (first CENTER, non-heading paragraph) -> '# '
  - Heading 1/2/3 -> '## ' / '### ' / '#### '
  - List Bullet -> '- '
  - EndNote Bibliography -> plain reference lines (under '## References')
  - Inline **bold** / *italic* / ^superscript^ (numbered citations)
  - Tables -> markdown tables (first row = header)
  - Embedded figures -> '![Figure N](figures/figureN.png)' using the build map
"""
import re
from docx import Document
from docx.oxml.ns import qn
from docx.text.paragraph import Paragraph
from docx.table import Table
from docx.enum.text import WD_ALIGN_PARAGRAPH

SRC = r"C:/D/workbuddy/科研/虚拟敲除/gate1/rewiring_study/manuscript_v7.docx"
OUT = r"C:/D/workbuddy/科研/虚拟敲除/gate1/rewiring_study/manuscript_v7.md"

FIG_MAP = {
    "1": "figure1_tsc_gnn_conceptual_framework_v4.png",
    "2": "figure2_prediction_benchmark.png",
    "3": "figure3_rewiring_heatmaps.png",
    "4": "figure4_cross_species.png",
    "5": "figure5_drug_reversal.png",
    "6": "figure6_L5_triangulation.png",
    "7": "figure7_evidence_ladder.png",
    "8": "cellchat_rewiring_heatmap.png",
    "S1": "figureS1_pc_correction.png",
}

MD_LVL = {1: 2, 2: 3, 3: 4}  # docx Heading N -> markdown '#'*k


def inline_md(runs):
    out = []
    cur = {"b": False, "i": False, "s": False}

    def close_unwanted(d):
        s = ""
        if cur["s"] and not d["s"]:
            s += "^"; cur["s"] = False
        if cur["i"] and not d["i"]:
            s += "*"; cur["i"] = False
        if cur["b"] and not d["b"]:
            s += "**"; cur["b"] = False
        return s

    def open_desired(d):
        s = ""
        if d["b"] and not cur["b"]:
            s += "**"; cur["b"] = True
        if d["i"] and not cur["i"]:
            s += "*"; cur["i"] = True
        if d["s"] and not cur["s"]:
            s += "^"; cur["s"] = True
        return s

    for r in runs:
        txt = r.text
        if not txt:
            continue
        d = {"b": bool(r.bold), "i": bool(r.italic), "s": bool(r.font.superscript)}
        out.append(close_unwanted(d) + open_desired(d) + txt)
    tail = ""
    if cur["s"]:
        tail += "^"; cur["s"] = False
    if cur["i"]:
        tail += "*"; cur["i"] = False
    if cur["b"]:
        tail += "**"; cur["b"] = False
    return "".join(out) + tail


def para_to_md(p, is_title=False):
    style = (p.style.name if p.style else "") or ""
    if style.startswith("Heading"):
        lvl = int(re.sub(r"\D", "", style) or 1)
        return "#" * MD_LVL.get(lvl, 2) + " " + inline_md(p.runs)
    if style == "List Bullet":
        return "- " + inline_md(p.runs)
    if style == "EndNote Bibliography":
        return p.text.replace("\t", " ").strip()
    if is_title:
        return "# " + inline_md(p.runs)
    return inline_md(p.runs)


def table_to_md(t):
    rows = []
    for row in t.rows:
        cells = []
        for c in row.cells:
            txt = " ".join(inline_md(pp.runs) for pp in c.paragraphs if pp.text.strip())
            txt = txt.replace("|", "/")
            cells.append(txt)
        rows.append(cells)
    if not rows:
        return ""
    ncols = len(rows[0])
    lines = ["| " + " | ".join(rows[0]) + " |",
             "| " + " | ".join(["---"] * ncols) + " |"]
    for r in rows[1:]:
        lines.append("| " + " | ".join(r) + " |")
    return "\n".join(lines)


def figure_ref_for(child, doc):
    nxt = child.getnext()
    if nxt is not None and nxt.tag == qn("w:p"):
        np = Paragraph(nxt, doc)
        m = re.search(r"Figure\s+(\d+|S\d+)", np.text or "")
        if m:
            num = m.group(1)
            fn = FIG_MAP.get(num, f"figure{num}.png")
            return f"![Figure {num}](figures/{fn})"
    return "<!-- [embedded figure] -->"


def main():
    doc = Document(SRC)
    body = doc.element.body
    lines = []
    title_emitted = False
    refs_heading_done = False
    prev_blank = False
    for child in body.iterchildren():
        if child.tag == qn("w:p"):
            p = Paragraph(child, doc)
            has_drawing = child.find(".//" + qn("w:drawing")) is not None
            if has_drawing:
                lines.append(figure_ref_for(child, doc))
                prev_blank = False
                continue
            style = (p.style.name if p.style else "") or ""
            if style == "EndNote Bibliography" and not refs_heading_done:
                lines.append("## References")
                refs_heading_done = True
            is_title = (not title_emitted) and p.alignment == WD_ALIGN_PARAGRAPH.CENTER \
                and not (p.style and p.style.name.startswith("Heading"))
            if is_title:
                title_emitted = True
            md = para_to_md(p, is_title=is_title)
            if md == "":
                if not prev_blank:
                    lines.append("")
                    prev_blank = True
            else:
                lines.append(md)
                prev_blank = False
        elif child.tag == qn("w:tbl"):
            t = Table(child, doc)
            md = table_to_md(t)
            if md:
                lines.append("")
                lines.append(md)
                lines.append("")
                prev_blank = True
    # strip leading/trailing blank lines
    while lines and lines[0] == "":
        lines.pop(0)
    while lines and lines[-1] == "":
        lines.pop()
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("Wrote", OUT)
    print("total lines:", len(lines))


if __name__ == "__main__":
    main()
