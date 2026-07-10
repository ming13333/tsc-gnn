"""分步计时诊断：定位 load_integrated_timeseries 的卡点。"""
import time, os
import numpy as np
import scanpy as sc
from gate1 import preprocessing as P
from gate1.data_acquisition import parse_series_matrix

t0 = time.time()
def tic(msg):
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)

# ---- 1) parse_series_matrix（网络）----
tic(">>> parse_series_matrix(GSE174574)")
samples = parse_series_matrix("GSE174574")
tic(f"<<< parse done: {len(samples)} samples")

# ---- 2) 逐个读 GSE174574 mtx（mmread 计时）----
ads = []
for s in samples:
    gsm_dir = os.path.join("data", "GSE174574", s["gsm"])
    paths = P._find_sample_paths(gsm_dir)
    if paths is None:
        tic(f"  skip {s['gsm']} (no mtx)"); continue
    ta = time.time()
    ad = P._read_10x(paths)
    tic(f"  read {s['gsm']} ({s['condition']}) cells={ad.n_obs} genes={ad.n_vars} in {time.time()-ta:.1f}s")
    ad.obs["sample"] = s["gsm"]
    ad.obs["condition"] = s["condition"]
    ad.obs["time_label"] = "0" if s["condition"] == "sham" else "1"
    ads.append(ad)

# ---- 3) concat + qc + normalize（load_gse174574_raw 后半）----
tic(">>> GSE174574 concat + qc + normalize")
for ad in ads:
    ad.var_names_make_unique()
a174 = sc.concat(ads, join="outer")
a174.var_names_make_unique(); a174.obs_names_make_unique()
a174.var["mt"] = a174.var_names.str.startswith("mt-") | a174.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(a174, qc_vars=["mt"], inplace=True)
a174 = a174[a174.obs["n_genes_by_counts"] > 200].copy()
sc.pp.filter_genes(a174, min_cells=3)
sc.pp.normalize_total(a174, target_sum=1e4)
sc.pp.log1p(a174)
tic(f"<<< GSE174574 ready: cells={a174.n_obs} genes={a174.n_vars}")

# ---- 4) load_gse225948（计时）----
tic(">>> load_gse225948(brain, male, W8, MCAO)")
a225 = P.load_gse225948("data", tissue="brain", sex="male", age="W8", condition="MCAO")
tic(f"<<< GSE225948 ready: cells={a225.n_obs} genes={a225.n_vars}")
tic("    GSE225948 time_label: " + str(a225.obs["time_label"].value_counts().to_dict()))

# ---- 5) align_genes ----
tic(">>> align_genes")
a174_a, a225_a, common = P.align_genes(a174, a225)
tic(f"<<< align done: common genes={len(common)}")

# ---- 6) concat 两集 ----
tic(">>> sc.concat([a174_a, a225_a])")
adata = sc.concat([a174_a, a225_a], join="outer")
adata.var_names_make_unique(); adata.obs_names_make_unique()
adata.obs["time_label"] = adata.obs["time_label"].astype(str)
tic(f"<<< integrated: cells={adata.n_obs} genes={adata.n_vars}")
tic("    study dist: " + str(adata.obs["study"].value_counts().to_dict()))
tic("    time dist: " + str(adata.obs["time_label"].value_counts().to_dict()))
print("DIAG_DONE", flush=True)
