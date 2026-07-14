# cellchat_lrdb.py
# 精选的 CellChatDB-mouse 配体-受体对子集（分泌信号 / 表面接触 / ECM-受体），
# 优先收录 CNS-卒中相关通路（小胶-神经-星形-内皮-外周免疫）。
# 仅用于按 CellChat 打分逻辑在 Python 中重实现细胞通讯分析（本机无 R）。
# 每条: (ligand, receptor, pathway)

LR_DB = [
    # ---- 趋化因子 / 细胞因子（分泌）----
    ("Ccl2", "Ccr2", "Chemokine"),
    ("Ccl3", "Ccr1", "Chemokine"),
    ("Ccl3", "Ccr5", "Chemokine"),
    ("Ccl4", "Ccr5", "Chemokine"),
    ("Ccl5", "Ccr5", "Chemokine"),
    ("Ccl7", "Ccr2", "Chemokine"),
    ("Ccl8", "Ccr2", "Chemokine"),
    ("Ccl20", "Ccr6", "Chemokine"),
    ("Ccl22", "Ccr4", "Chemokine"),
    ("Cxcl1", "Cxcr2", "Chemokine"),
    ("Cxcl2", "Cxcr2", "Chemokine"),
    ("Cxcl9", "Cxcr3", "Chemokine"),
    ("Cxcl10", "Cxcr3", "Chemokine"),
    ("Cxcl12", "Cxcr4", "Chemokine"),
    ("Cxcl16", "Cxcr6", "Chemokine"),
    ("Cx3cl1", "Cx3cr1", "Chemokine"),
    ("Csf1", "Csf1r", "Cytokine"),
    ("Csf2", "Csf2ra", "Cytokine"),
    ("Il34", "Csf1r", "Cytokine"),
    ("Il1b", "Il1r1", "Cytokine"),
    ("Il6", "Il6ra", "Cytokine"),
    ("Il10", "Il10ra", "Cytokine"),
    ("Il4", "Il4ra", "Cytokine"),
    ("Lif", "Lifr", "Cytokine"),
    ("Osm", "Osmr", "Cytokine"),
    ("Tnf", "Tnfrsf1a", "Cytokine"),
    ("Tgfb1", "Tgfbr1", "TGFb"),
    ("Tgfb1", "Tgfbr2", "TGFb"),
    ("Tgfb2", "Tgfbr2", "TGFb"),
    ("Tgfb3", "Tgfbr2", "TGFb"),
    ("Ifng", "Ifngr1", "Cytokine"),
    ("Vegfa", "Kdr", "VEGF"),
    ("Vegfa", "Flt1", "VEGF"),
    ("Vegfc", "Flt4", "VEGF"),
    # ---- 生长因子 / 形态发生原（分泌）----
    ("Pdgfa", "Pdgfra", "PDGF"),
    ("Pdgfb", "Pdgfra", "PDGF"),
    ("Egfr", "Egfr", "EGF"),
    ("Hgf", "Met", "HGF"),
    ("Bmp4", "Bmpr1a", "BMP"),
    ("Bmp7", "Bmpr1a", "BMP"),
    ("Tgfb1", "Acvrl1", "TGFb"),
    ("Wnt5a", "Fzd1", "Wnt"),
    ("Wnt3a", "Fzd1", "Wnt"),
    ("Wnt7a", "Fzd5", "Wnt"),
    ("Spp1", "Cd44", "Osteopontin"),
    ("Spp1", "Itgav", "Osteopontin"),
    ("Gas6", "Axl", "Gas6"),
    ("Thbs1", "Cd47", "Thrombospondin"),
    ("Thbs1", "Cd36", "Thrombospondin"),
    ("Thbs2", "Cd47", "Thrombospondin"),
    # ---- 补体 / 其他分泌 ----
    ("C1qa", "C1qb", "Complement"),
    ("C3", "C3ar1", "Complement"),
    ("C5", "C5ar1", "Complement"),
    # ---- 细胞-细胞接触 / 表面（免疫突触、黏附）----
    ("Icam1", "Itgal", "Adhesion"),
    ("Icam1", "Itgb2", "Adhesion"),
    ("Vcam1", "Itga4", "Adhesion"),
    ("Vcam1", "Itgav", "Adhesion"),
    ("Selp", "Selplg", "Adhesion"),
    ("Sele", "Selplg", "Adhesion"),
    ("Cd40", "Cd40lg", "Costim"),
    ("Cd70", "Cd27", "Costim"),
    ("Icoslg", "Icos", "Costim"),
    ("Pecam1", "Pecam1", "Adhesion"),
    ("Ncam1", "Ncam1", "Adhesion"),
    ("Cdh1", "Cdh1", "Adhesion"),
    ("Cdh2", "Cdh2", "Adhesion"),
    ("Efnb2", "Ephb2", "Ephrin"),
    ("Efnb3", "Ephb3", "Ephrin"),
    ("Epha2", "Efna1", "Ephrin"),
    ("Jaml", "Jam3", "Adhesion"),
    ("Tigit", "Pvr", "Checkpoint"),
    ("Havcr2", "Havcr1", "Checkpoint"),
    ("Pdcd1", "Cd274", "Checkpoint"),
    ("Ctla4", "Cd80", "Checkpoint"),
    # ---- 存活 / 清除 ----
    ("Mif", "Cd74", "CD74"),
    ("App", "Sorl1", "Amyloid"),
    ("S100a9", "Tlr4", "Alarmin"),
    ("S100a8", "Tlr4", "Alarmin"),
    ("Hmgb1", "Tlr4", "Alarmin"),
    # ---- 神经相关递质/受体 ----
    ("Gdnf", "Ret", "GDNF"),
    ("Bdnf", "Ntrk2", "Neurotrophin"),
    ("Ntf3", "Ntrk3", "Neurotrophin"),
    ("Nrg1", "Erbb4", "Neuregulin"),
    ("Sema3a", "Nrp1", "Semaphorin"),
    ("Sema4d", "Plexdc1", "Semaphorin"),
    ("Slit2", "Robo1", "Slit"),
    # ---- ECM-受体 ----
    ("Col1a1", "Cd44", "ECM"),
    ("Col4a1", "Cd44", "ECM"),
    ("Fn1", "Itga5", "ECM"),
    ("Fn1", "Itgb1", "ECM"),
    ("Lama2", "Itga7", "ECM"),
    ("Lamb2", "Itgb1", "ECM"),
    ("Tnc", "Itgav", "ECM"),
    ("Vcan", "Cd44", "ECM"),
    ("Sdc1", "Cxcl12", "ECM"),
    ("Hpse", "Sdc1", "ECM"),
    # ---- 代谢/其他通路 ----
    ("Apoe", "Ldlr", "Lipoprotein"),
    ("Apoe", "Lrp1", "Lipoprotein"),
    ("Cxcl10", "Sdc1", "Chemokine"),
    ("Pf4", "Pf4r", "Chemokine"),
    ("Ccl2", "Sdc1", "Chemokine"),
    ("Tgfa", "Egfr", "EGF"),
    ("Areg", "Egfr", "EGF"),
    ("Fgf1", "Fgfr1", "FGF"),
    ("Fgf2", "Fgfr1", "FGF"),
    ("Igf1", "Igf1r", "IGF"),
    ("Inhba", "Acvr1b", "TGFb"),
    ("Pdgfa", "Pdgfrb", "PDGF"),
    ("Tgfb2", "Acvrl1", "TGFb"),
]


def lr_pairs_unique():
    """返回去重的 (ligand, receptor) 集合。"""
    return sorted({(l, r) for l, r, _ in LR_DB})


def pathway_of(ligand, receptor):
    for l, r, p in LR_DB:
        if l == ligand and r == receptor:
            return p
    return "Unknown"
