#!/usr/bin/env python3
"""Human adult temporal-lobe surgical-control brain endothelial RNA reference.

Original GEO GSE256490, GSM8101534--GSM8101539 (TL_patient1_EC .. 6_EC).
Six separately deposited CD31+CD45- FACS endothelial preparations, NOT six
ischemic stroke subjects; control tissue includes epilepsy surgical setting.
No differential expression, cell-level p values, PET target validation or
unverified single-cell author annotation transferred.
"""
from __future__ import annotations
import csv,gzip,hashlib,json,math,re,tempfile,time,urllib.request
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"human_vascular_reference"
SOURCE=OUT/"19_GSE256490_adult_TL_FACS_EC_candidate_libraries.csv"
PANEL=ROOT/"config"/"candidate_catalog.csv"
EC_MARKERS=["PECAM1","CLDN5","KDR","VWF","CDH5","EMCN","ESAM"]
OTHER_MARKERS=["PTPRC","AQP4","PDGFRB","ACTA2","P2RY12","MOBP"]
MAX_COMPRESSED_FILE=180*1024*1024

def require(expr,message):
    if not expr:raise AssertionError(message)

def table(path):
    with path.open("r",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def write(name,rows):
    require(bool(rows),"Empty report "+name)
    with (OUT/name).open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="ignore")
        writer.writeheader();writer.writerows(rows)

def get(url,path):
    require(url.startswith("https://ftp.ncbi.nlm.nih.gov/geo/samples/"),
            "Do not download from non-source endpoint")
    require(url.endswith((".mtx.gz",".tsv.gz")),"Unexpected GEO data type")
    last=None
    for attempt in range(2):
        try:
            q=urllib.request.Request(url,headers={"User-Agent":"IMACS-vascular-reference/1.0"})
            with urllib.request.urlopen(q,timeout=180) as r,path.open("wb") as f:
                length=r.headers.get("content-length")
                if length is not None:require(int(length)<=MAX_COMPRESSED_FILE,
                                              "Unexpected oversized source file")
                total=0
                while chunk:=r.read(1024*1024):
                    total+=len(chunk)
                    require(total<=MAX_COMPRESSED_FILE,"Source file exceeded safe download bound")
                    f.write(chunk)
            require(path.stat().st_size>100,"Empty source download")
            return hashlib.sha256(path.read_bytes()).hexdigest(),path.stat().st_size
        except (OSError,TimeoutError) as exc:
            last=exc
            path.unlink(missing_ok=True)
            if attempt==0:time.sleep(3)
    raise RuntimeError("GEO source download blocked: "+str(last))

def extract_sample(gsm,record,genes,outrows,markerrows,quality,checks,folder):
    data={}
    for kind in ("matrix","features","barcodes"):
        url=record[f"{kind}_url"]
        require(url and "/"+gsm+"/suppl/"+gsm+"_" in url,
                f"Unexpected source GEO file for {gsm} {kind}")
        p=folder/(gsm+"_"+kind+".gz")
        sha,size=get(url,p)
        data[kind]=p
        checks.append(dict(gsm=gsm,title=record["title"],kind=kind,
                           source_url=url,bytes=size,sha256=sha))
    with gzip.open(data["features"],"rt",encoding="utf-8") as f:
        feat=pd.read_csv(f,sep="\t",header=None,dtype=str)
    require(2<=feat.shape[1]<=4,f"Unexpected feature columns {gsm}")
    symbols=feat.iloc[:,1].fillna("").str.upper().to_numpy()
    require(len(set(feat.iloc[:,0]))==len(feat),"Duplicated Ensembl gene feature ID")
    with gzip.open(data["barcodes"],"rt") as f:
        barcodes=[s.strip() for s in f if s.strip()]
    require(len(set(barcodes))==len(barcodes),"Duplicate barcodes")
    with gzip.open(data["matrix"],"rb") as f:
        matrix=mmread(f).tocsr()
    require(matrix.shape==(len(symbols),len(barcodes)),
            f"Matrix dimensions differ from exact GEO sources {gsm}")
    require(sparse.issparse(matrix),"Unexpected dense library")
    require(np.isfinite(matrix.data).all() and (matrix.data>=0).all(),
            "Invalid negative or nonfinite expression")
    # If not a UMI-like integer source, record it. We use only descriptive RNA.
    fraction_nonint=float(np.mean(np.abs(matrix.data-np.rint(matrix.data))>1e-8)) if matrix.nnz else 0.
    X=matrix.transpose().tocsr()
    total=np.asarray(X.sum(axis=1)).ravel()
    ngenes=np.diff(X.indptr)
    mt_cols=np.flatnonzero(np.char.startswith(symbols.astype(str),"MT-"))
    mt=np.asarray(X[:,mt_cols].sum(axis=1)).ravel() if len(mt_cols) else np.zeros(X.shape[0])
    pct_mt=100*mt/np.maximum(total,1)
    # Track original cell capture, apply *per preparation* conservative thresholds.
    accept=((ngenes>=200)&(pct_mt<20)&
            (ngenes<=np.quantile(ngenes,.995))&
            (total<=np.quantile(total,.995)))
    require(accept.sum()>=50,f"Insufficient captured nuclei/cells {gsm}")
    X=X[accept].tocsr()
    total=total[accept]
    quality.append(dict(gsm=gsm,title=record["title"],
         source="GSE256490",
         source_includes_epilepsy_surgical_controls=bool("Epilepsy" in record["matrix_url"]),
         captured_cells_before=int(len(barcodes)),qc_accepted_cells=int(X.shape[0]),
         qc_retention_fraction=float(accept.mean()),
         median_features_after=float(np.median(ngenes[accept])),
         median_UMI_after=float(np.median(total)),
         median_mito_fraction_after=float(np.median(pct_mt[accept])),
         noninteger_original_matrix_fraction=fraction_nonint,
         measurement_level="FACS CD31+ CD45- pooled-cell distribution, not donor-level stroke comparison"))
    idx=defaultdict(list)
    chosen=set(genes+EC_MARKERS+OTHER_MARKERS)
    for j,s in enumerate(symbols):
        if s in chosen:idx[s].append(j)
    for gene in sorted(chosen):
        columns=idx.get(gene,[])
        if not columns:
            mean=""; detected=""; median=""
        else:
            vals=np.asarray(X[:,columns].sum(axis=1)).ravel()
            normalized=np.log1p(vals/np.maximum(total,1)*1e4)
            mean=float(normalized.mean())
            detected=float(np.mean(vals>0))
            median=float(np.median(normalized))
        row=dict(gsm=gsm,source_title=record["title"],gene=gene,
                 n_original_symbol_features=len(columns),captured_cells_after_QC=X.shape[0],
                 fraction_cells_detected=detected,
                 average_cell_log1p_10000=mean,
                 median_cell_log1p_10000=median,
                 original_symbol_present=bool(columns),
                 readout="DESCRIPTIVE_SOURCE_GEO_FACS_EC_SURGICAL_CONTROL")
        if gene in genes:outrows.append(row)
        else:markerrows.append(row)
    print("HUMAN_EC_SAMPLE",gsm,record["title"],"cells",len(barcodes),"->",X.shape[0],
          "genes",len(symbols),"nonint",fraction_nonint,flush=True)
    del matrix,X

def main():
    source=table(SOURCE)
    require(len(source)==6,"Expect six named adult temporal-lobe EC patients")
    expected={f"TL_patient{i}_EC" for i in range(1,7)}
    require({r["title"] for r in source}==expected,"Human source identifiers changed")
    require(len(set(r["gsm"] for r in source))==6,"Missing GSM library")
    require(all(r["candidate_adult_temporal_lobe_FACS_EC"]=="True" and
                "treatment: control" in r["characteristics"] and
                "tissue: adult/control brain" in r["characteristics"]
                for r in source),"Noncontrol/unsupported human brain source")
    catalog=table(PANEL)
    genes=list(dict.fromkeys(r["gene"].upper() for r in catalog))
    require(len(genes)==29,"Panel gene count changed")
    rows=[];markers=[];quality=[];checks=[]
    with tempfile.TemporaryDirectory() as d:
        directory=Path(d)
        for sample in source:
            extract_sample(sample["gsm"],sample,genes,rows,markers,quality,checks,directory)
            for p in directory.iterdir(): p.unlink()
    require(len(rows)==6*29,"Incomplete six-source 29-target expression matrix")
    summarised=[]
    for gene in genes:
        gr=[r for r in rows if r["gene"]==gene]
        mustread=[r for r in gr if r["original_symbol_present"]]
        fractions=[r["fraction_cells_detected"] for r in mustread]
        summarised.append(dict(
            gene=gene,n_control_preparations=6,
            n_source_symbol_covered=len(mustread),
            n_preparations_with_detection_ge_10pct=
               sum(x>=0.10 for x in fractions),
            median_fraction_cells_detected_across_covered_sources=
               float(np.median(fractions)) if fractions else "",
            range_fraction_cells_detected_across_covered_sources=
               f"{min(fractions):.4f}--{max(fractions):.4f}" if fractions else "",
            human_baseline_interpretation=
               "POSTMORTEM_BRAIN_STROKE_EXPRESSION_NOT_MEASURED_HERE",
            vascular_selectivity="NO_UNSORTED_COMPARATOR_IN_THIS_ANALYSIS",
            PET_targetability="NOT_INFERABLE_FROM_RNA_DETECTION"))
    OUT.mkdir(parents=True,exist_ok=True)
    write("20_GSE256490_human_TL_surgical_control_FACS_EC_expression_by_preparation.csv",rows)
    write("20_GSE256490_human_TL_candidate_detection_across_six_preparations.csv",summarised)
    write("20_GSE256490_human_TL_lineage_marker_expression_by_preparation.csv",markers)
    write("20_GSE256490_human_TL_original_matrix_QC.csv",quality)
    write("20_GSE256490_human_TL_source_sha256.csv",checks)
    manifest={
      "accession":"GSE256490","source_publication":"Wälchli et al. Nature 2024",
      "doi":"10.1038/s41586-024-07493-y",
      "specimen":"six distinct GEO TL_patient[1-6]_EC surgical adult temporal lobe samples",
      "sampling":"FACS sorted CD31-positive CD45-negative endothelial-enriched tissue; may contain other cell types",
      "donor_context":"controls from brain surgery; first two matrix filenames explicitly say Epilepsy; not healthy population controls",
      "n_original_GEO_preparations":len(source),
      "GSM":[r["gsm"] for r in source],
      "individual_cells_biological_replicates":False,
      "normalization":"per-cell library-size normalized log1p(10000*source_U MI/total); expression and detection descriptive",
      "cross_stroke_comparison":"NOT PERFORMED: NO ISCHEMIC BRAIN IN THIS SELECTION",
      "interpretation":"Patient-level baseline detection in human surgical brain EC, NOT independent validation of stroke induction or human BBB luminal protein",
      "gene_symbol_mapping":"exact same deposited HGNC symbol (case insensitive), duplicate source features summed; absent != biologically absent",
      "source_sha256_in":"20_GSE256490_human_TL_source_sha256.csv",
      "protein_or_PET_feasibility":False,
      "summary":{
        "n_target_genes":len(genes),
        "genes_covered_in_all_six_preparations":sum(r["n_source_symbol_covered"]==6 for r in summarised),
        "genes_detected_in_at_least_10pct_cells_in_all_six":sum(r["n_preparations_with_detection_ge_10pct"]==6 for r in summarised),
        "total_cells_after_QC":sum(r["qc_accepted_cells"] for r in quality),
      }
    }
    (OUT/"20_human_TL_brain_endothelial_baseline_manifest.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print("HUMAN_EC_BASELINE PASS",json.dumps(manifest["summary"]),flush=True)

if __name__=="__main__":main()
