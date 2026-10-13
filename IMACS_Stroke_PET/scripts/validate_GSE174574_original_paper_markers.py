#!/usr/bin/env python3
"""Replicate original IMACS GSE174574 classifier and assess ORIGINAL-PAPER markers.
NO author cell-ID labels supplied; this tests cross-marker concordance, not exact author reannotation.
"""

import os, sys, gzip, hashlib, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import requests

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=90, facecolor="white")

ROOT = Path(os.environ.get("IMACS_ROOT", "/content/imacs_stroke_pet" if Path("/content").exists() else "imacs_stroke_pet"))
RAW = ROOT/"raw"
RESULTS = ROOT/"results"
FIGS = ROOT/"figures"
for p in (RAW, RESULTS, FIGS): p.mkdir(parents=True, exist_ok=True)

print("ROOT:", ROOT)
print("Python:", sys.version.split()[0], "scanpy:", sc.__version__)

SAMPLES = pd.DataFrame([
    ["GSM5319987","sham1","Sham"],
    ["GSM5319988","sham2","Sham"],
    ["GSM5319989","sham3","Sham"],
    ["GSM5319990","MCAO1","MCAO"],
    ["GSM5319991","MCAO2","MCAO"],
    ["GSM5319992","MCAO3","MCAO"],
], columns=["gsm","prefix","condition"])
print(SAMPLES.to_string(index=False))

def geo_block(gsm):
    return gsm[:7] + "nnn"

def geo_url(gsm, prefix, filename):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{geo_block(gsm)}/{gsm}/suppl/{gsm}_{prefix}_{filename}"

def download(url, path):
    path=Path(path)
    if path.exists() and path.stat().st_size>0: return path
    with requests.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        with open(path,"wb") as f:
            for chunk in r.iter_content(1<<20):
                if chunk: f.write(chunk)
    return path

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1<<20), b""): h.update(chunk)
    return h.hexdigest()

def fetch_sample(gsm,prefix,sample_id,condition):
    d=RAW/gsm; d.mkdir(exist_ok=True)
    paths={}
    for fn in ["matrix.mtx.gz","barcodes.tsv.gz","genes.tsv.gz"]:
        u=geo_url(gsm,prefix,fn)
        p=download(u,d/f"{gsm}_{prefix}_{fn}")
        paths[fn]=p
    with gzip.open(paths["matrix.mtx.gz"],"rb") as f:
        X=mmread(f).tocsr().T
    with gzip.open(paths["barcodes.tsv.gz"],"rt") as f:
        barcodes=[x.rstrip() for x in f]
    feats=pd.read_csv(paths["genes.tsv.gz"],sep="\t",header=None,compression="gzip")
    symbols=feats.iloc[:,1].astype(str).values if feats.shape[1]>1 else feats.iloc[:,0].astype(str).values
    a=ad.AnnData(X=X)
    a.obs_names=[f"{sample_id}:{b}" for b in barcodes]
    a.var_names=pd.Index(symbols)
    a.var_names_make_unique()
    a.obs["sample_id"]=sample_id
    a.obs["condition"]=condition
    return a, paths

objs=[]; prov=[]
for r in SAMPLES.itertuples(index=False):
    print("Loading",r.gsm,r.prefix)
    a,paths=fetch_sample(r.gsm,r.prefix,r.prefix,r.condition)
    objs.append(a)
    for kind,p in paths.items():
        prov.append({"gsm":r.gsm,"sample_id":r.prefix,"condition":r.condition,"file":p.name,"bytes":p.stat().st_size,"sha256":sha256(p)})

pd.DataFrame(prov).to_csv(RESULTS/"01_input_checksums.csv",index=False)
adata=ad.concat(objs,join="outer",merge="same",index_unique=None)
assert adata.obs_names.is_unique, "Duplicate cell identifiers across samples!"
print(adata)

adata.var["mt"]=adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata,qc_vars=["mt"],inplace=True)

before=adata.obs.groupby("sample_id").agg(
    n_cells=("sample_id","size"),
    median_genes=("n_genes_by_counts","median"),
    median_umis=("total_counts","median"),
    median_pct_mt=("pct_counts_mt","median"))


keep=(adata.obs["n_genes_by_counts"].values>=200) & (adata.obs["pct_counts_mt"].values<20)
for s,names in adata.obs.groupby("sample_id").groups.items():
    idx=adata.obs_names.get_indexer(names)
    d=adata.obs.iloc[idx]
    keep[idx] &= d["n_genes_by_counts"].values <= d["n_genes_by_counts"].quantile(.995)
    keep[idx] &= d["total_counts"].values <= d["total_counts"].quantile(.995)

adata=adata[keep].copy()
adata.layers["counts"]=adata.X.copy()
after=adata.obs.groupby("sample_id").size().rename("n_cells_after")
audit=before.join(after)

audit.to_csv(RESULTS/"02_QC_audit.csv")

MARKERS={
"Endothelial":["Cldn5","Pecam1","Kdr","Flt1","Esam","Emcn","Ramp2","Klf2"],
"Pericyte":["Pdgfrb","Rgs5","Kcnj8","Abcc9","Anpep","Cd248","Cspg4"],
"VSMC":["Acta2","Tagln","Myh11","Cnn1","Rgs5"],
"Astrocyte":["Aqp4","Slc1a3","Aldoc","Gja1","S100b","Aldh1l1"],
"Microglia":["P2ry12","Tmem119","C1qa","C1qb","Hexb","Cx3cr1"],
"Myeloid":["Ptprc","Lyz2","Tyrobp","Aif1","Fcgr3"],
"Neuron":["Rbfox3","Snap25","Syt1","Map2","Tubb3"],
"Oligodendrocyte":["Mbp","Plp1","Mog","Mag","Cnp"],
"OPC":["Pdgfra","Olig1","Olig2","Sox10"],
"VLMC_Fibroblast":["Col1a1","Col1a2","Dcn","Lum","Col3a1"],
"Ependymal_CP":["Ttr","Foxj1","Clic6"]
}

sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
adata.raw=adata.copy()  # Restore original pipeline object state before clustering

sc.pp.highly_variable_genes(adata,n_top_genes=2500,flavor="seurat",batch_key="sample_id")
sc.tl.pca(adata,n_comps=40,use_highly_variable=True)
sc.pp.neighbors(adata,n_neighbors=15,n_pcs=30)
sc.tl.umap(adata)  # Preserve original pipeline's computation order before Leiden
sc.tl.leiden(adata,resolution=.6,key_added="leiden")

for lineage,genes in MARKERS.items():
    genes=[g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata,genes,score_name=f"score_{lineage}",use_raw=True)

score_cols=[f"score_{x}" for x in MARKERS]
cs=adata.obs.groupby("leiden")[score_cols].mean()
cs["best"]=cs.idxmax(axis=1).str.replace("score_","",regex=False)
vals=np.sort(cs[score_cols].values,axis=1)
cs["margin"]=vals[:,-1]-vals[:,-2]
cs["broad_celltype"]=np.where(cs["margin"]>=.05,cs["best"],"Ambiguous")
adata.obs["broad_celltype"]=adata.obs["leiden"].map(cs["broad_celltype"]).astype("category")


cs.to_csv(RESULTS/"03_cluster_annotation_scores.csv")



# Original pipeline's source QC, UMAP, Leiden, score_genes and broad calls are reproduced in the same order.\n# Source-independent marker cross-check from Zheng et al., JCBFM 2022
# DOI 10.1177/0271678X211026770: Itm2a/EC, Hexb/MG, Aldoc/astrocyte,
# Plp1/oligodendrocyte, Ttr/ependymal, Acta2/VSMC.
# This is NOT recovery of original per-cell authors' annotations.
import json
from collections import defaultdict
REF_MARKERS={
    "Itm2a":"Endothelial",
    "Hexb":"Microglia",
    "Aldoc":"Astrocyte",
    "Plp1":"Oligodendrocyte",
    "Ttr":"Ependymal_CP",
    "Acta2":"VSMC"
}
assert list(SAMPLES["gsm"])==[
    "GSM5319987","GSM5319988","GSM5319989","GSM5319990","GSM5319991","GSM5319992"
]
expected_prov=pd.read_csv(
    "IMACS_Stroke_PET/results/primary/01_input_checksums.csv"
)
given={(r.gsm,r.file):r.sha256 for r in expected_prov.itertuples(index=False)}
observed={(r["gsm"],r["file"]):r["sha256"] for r in prov}
assert given==observed, "Source checksums mismatch: original marker audit MUST stop"
existing=pd.read_csv(
    "IMACS_Stroke_PET/results/primary/04_celltype_counts_by_sample.csv"
).set_index("sample_id")
new_counts=pd.crosstab(adata.obs["sample_id"],adata.obs["broad_celltype"])
old_aligned=existing.reindex(index=SAMPLES.prefix,columns=existing.columns,fill_value=0)
new_aligned=new_counts.reindex(index=SAMPLES.prefix,columns=existing.columns,fill_value=0)
class_drift=(new_aligned-old_aligned).abs()
n_disagree=int((class_drift>0).sum().sum())
max_abs=int(class_drift.to_numpy().max())
print("CLUSTER_COUNT_DRIFT",n_disagree,"changed sample-class slots; max abs diff",max_abs)
# Avoid treating marker identity as full, original author-state reconstruction.
counts=adata.layers["counts"]
def sparse_vector(matrix, idx, j):
    val=matrix[idx,j]
    return np.asarray(val.toarray()).ravel() if sp.issparse(val) else np.asarray(val).ravel()
rows=[]
for sample in SAMPLES.prefix:
    for lineage in MARKERS:
        idx=np.flatnonzero(
            (adata.obs["sample_id"].to_numpy()==sample) &
            (adata.obs["broad_celltype"].astype(str).to_numpy()==lineage)
        )
        for gene,target in REF_MARKERS.items():
            present=gene in adata.var_names
            if not len(idx):
                rows.append(dict(sample_id=sample,broad_celltype=lineage,
                    gene=gene,published_marker_expected_class=target,
                    n_cells=0,source_feature_present=present,
                    mean_log1p="",fraction_detected="",
                    marker_independent_of_classification_panel=(gene=="Itm2a")))
                continue
            if not present:
                rows.append(dict(sample_id=sample,broad_celltype=lineage,
                    gene=gene,published_marker_expected_class=target,
                    n_cells=len(idx),source_feature_present=False,
                    mean_log1p="",fraction_detected="",
                    marker_independent_of_classification_panel=(gene=="Itm2a")))
                continue
            j=adata.var_names.get_loc(gene)
            x=sparse_vector(adata.X,idx,j)
            c=sparse_vector(counts,idx,j)
            rows.append(dict(sample_id=sample,broad_celltype=lineage,
                gene=gene,published_marker_expected_class=target,
                n_cells=len(idx),source_feature_present=True,
                mean_log1p=float(x.mean()),
                fraction_detected=float(np.mean(c>0)),
                marker_independent_of_classification_panel=(gene=="Itm2a")))
marker_table=pd.DataFrame(rows)
marker_table.to_csv(
    RESULTS/"16_Zheng_original_paper_marker_concordance_by_sample.csv",index=False
)
marker_overview=[]
for marker,expected in REF_MARKERS.items():
    subset=marker_table[(marker_table.gene==marker)&(marker_table.n_cells>0)&
                        (marker_table.source_feature_present)]
    by_class=subset.groupby("broad_celltype").agg(
        n_samples=("sample_id","nunique"),
        average_sample_mean_log1p=("mean_log1p","mean"),
        average_sample_fraction_detected=("fraction_detected","mean"))
    by_class=by_class.sort_values("average_sample_mean_log1p",ascending=False)
    top=by_class.index.tolist()[:3]
    expected_n=int(by_class.loc[expected,"n_samples"]) if expected in by_class.index else 0
    marker_overview.append(dict(
        published_original_marker=marker,
        published_expected_class=expected,
        feature_present=marker in adata.var_names,
        expected_class_samples=expected_n,
        top_three_inhouse_broad_classes=";".join(top),
        expected_in_top_three=bool(expected in top),
        quantitative_source="mean of sample-wise normalized expression, descriptive",
        not_original_author_cell_id_reproduction=True
    ))
overview=pd.DataFrame(marker_overview)
overview.to_csv(
    RESULTS/"16_Zheng_author_marker_specificity_crosscheck.csv",index=False
)
report=dict(
    dataset="GSE174574",
    original_study_doi="10.1177/0271678X211026770",
    original_geo_samples=list(SAMPLES.gsm),
    biological_design="three deposited sham and three MCAO mouse sample labels",
    caveat="GEO rep identifiers do not expose a direct mapping from each cell to original raw animal beyond the six libraries",
    original_paper_annotations="17 original detailed types including capillary, arterial, and venous EC; no per-cell author label file was supplied in our source inputs",
    current_annotations="our 11 broad marker-score classes; not one-to-one comparable with published authors",
    independent_check_limit="Itm2a independent of original classifier marker list; other published markers overlap existing classifier panels",
    original_paper_marker_presence={k:bool(k in adata.var_names) for k in REF_MARKERS},
    expected_class_in_top_three={x["published_original_marker"]:bool(x["expected_in_top_three"]) for x in marker_overview},
    source_SHA256_matches_archived_six_samples=True,
    archived_compartment_counts_equal=bool(n_disagree==0),
    n_sample_by_celltype_count_slots_changed=n_disagree,
    maximal_absolute_change_in_class_counts=max_abs,
    marker_concordance="Descriptive only; no published original per-cell cluster transfer",
    inferential_outcome="UNCHANGED: no additional animal-level differential expression conducted",
    protein_PET="No transcript marker check proves protein targetability",
)
(RESULTS/"16_Zheng_original_marker_validation_manifest.json").write_text(json.dumps(
    report,indent=2)+"\n",encoding="utf-8")
print("MOUSE_MARKER_AUDIT PASS",json.dumps(report),flush=True)
