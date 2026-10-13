#!/usr/bin/env python3
"""Independent source-backed audit: GSE162955 six matched human stroke brains.

Uses genuine NCBI GEO processed series matrix (SHA256 pinned), archived
unambiguous GPL17586 probe-to-symbol map, human pair table and frozen panel.
Checks GEO GSM/title/age/sex pairing, original per-probe intensity, 64 exact
sign-flips, BH, and six leave-one-donor-out perturbations.
No cell-specific inference, no meta-analysis and no PET target ranking.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import itertools
import json
import math
import statistics
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results"
HUMAN=BASE/"human"
OUT=BASE/"human_replication_qc"
SOURCE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162955/matrix/GSE162955_series_matrix.txt.gz"
PATIENTS={
  "16":("GSM4968878","GSM4968879"),
  "22":("GSM4968880","GSM4968881"),
  "32":("GSM4968882","GSM4968883"),
  "33":("GSM4968884","GSM4968885"),
  "35":("GSM4968886","GSM4968887"),
  "36":("GSM4968888","GSM4968889"),
}
GSM=[s for pair in PATIENTS.values() for s in pair]
EXPECTED_SHA256="4cbb441eff1617bfac880804717dc8a2035c9e426b760e2eb0104ce6a2ce5497"

def must(ok,msg):
    if not ok: raise AssertionError(msg)

def read(path):
    with path.open("r",newline="",encoding="utf-8") as h:
        return list(csv.DictReader(h))

def save(path,rows):
    must(len(rows)>0,"No output for "+str(path))
    with path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def close(a,b,atol=1e-8):
    return math.isclose(float(a),float(b),rel_tol=atol,abs_tol=atol)

def meta_parse(path):
    meta=defaultdict(list)
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"): break
            if not line.startswith("!Sample_"): continue
            parts=next(csv.reader([line.rstrip("\r\n")],delimiter="\t",quotechar='"'))
            key=parts[0].split(" = ",1)[0]
            meta[key].append(parts[1:])
    for key in ("!Sample_geo_accession","!Sample_title"):
        must(key in meta, "Missing GEO metadata "+key)
    accession=meta["!Sample_geo_accession"][0]
    title=meta["!Sample_title"][0]
    must(len(accession)==12 and len(title)==12, "Expected 12 GEO patient samples")
    titles=dict(zip(accession,title))
    must(len(titles)==12 and set(titles)==set(GSM),"Unexpected GEO accessions or order")
    for patient,(ctrl,case) in PATIENTS.items():
        must(titles[ctrl].lower()==f"{patient}_contralateral" and
             titles[case].lower()==f"{patient}_infarct_core",
             f"Human GEO patient pair mismatch {patient} {titles[ctrl]} {titles[case]}")
    characteristics={}
    for row in meta.get("!Sample_characteristics_ch1",[]):
        for gsm,value in zip(accession,row):
            key, _, val=value.partition(": ")
            characteristics.setdefault(gsm,{})[key.strip().lower()]=val.strip()
    # Demographics should be consistent within a matched donor if present.
    for patient,(ctrl,case) in PATIENTS.items():
        for key in ("gender","age","sex"):
            x=characteristics.get(ctrl,{}).get(key)
            y=characteristics.get(case,{}).get(key)
            if x and y:
                must(x==y,f"Donor demographic mismatch {patient}: {key}")
    return titles,characteristics

def source_probe_extract(path,wanted):
    hits={}
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"): break
        header=next(csv.reader([next(f).strip()],delimiter="\t"))
        must(header[0]=="ID_REF" and set(header[1:])==set(GSM) and len(header)==13,
             "Series matrix sample headers mismatch")
        for line in f:
            if line.startswith("!series_matrix_table_end"):break
            values=next(csv.reader([line.strip()],delimiter="\t"))
            must(len(values)==13,"Damaged GEO matrix record")
            key=values[0]
            if key not in wanted:continue
            must(key not in hits,"Duplicate candidate probe in source matrix "+key)
            hits[key]={gsm:float(value) for gsm,value in zip(header[1:],values)}
    must(set(hits)==set(wanted),f"Missing {len(set(wanted)-set(hits))} target probes in GEO matrix")
    must(all(all(math.isfinite(v) for v in sample.values()) for sample in hits.values()),
         "Nonfinite GEO human candidate probe intensity")
    return hits

def signflip(deltas):
    must(len(deltas)==6 or len(deltas)==5,"Invalid pair count")
    observed=statistics.mean(deltas)
    null=[abs(statistics.mean(deltas[i]*signs[i] for i in range(len(deltas))))
          for signs in itertools.product((-1,1),repeat=len(deltas))]
    p=sum(x>=abs(observed)-1e-12 for x in null)/len(null)
    return observed,p

def bh(pvalues):
    indices=sorted(range(len(pvalues)),key=lambda i:pvalues[i])
    results=[None]*len(pvalues)
    last=1.0
    for k in range(len(indices),0,-1):
        i=indices[k-1]
        last=min(last,pvalues[i]*len(indices)/k)
        results[i]=last
    return results

def main():
    archived=read(HUMAN/"10_human_paired_core_vs_contralateral_candidate_screen.csv")
    archive_gene={x["gene"].upper():x for x in archived}
    must(len(archive_gene)==29 and len(archived)==29,"Frozen 29-gene panel mismatch")
    original_pairs=read(HUMAN/"10_human_per_patient_candidate_intensities.csv")
    original_pairs_by={(x["gene"].upper(),x["patient_id"]):x for x in original_pairs}
    must(len(original_pairs_by)==len(original_pairs),"Duplicate original gene/donor")
    # This map was extracted from GPL17586 by the earlier notebook.
    probe_rows=read(HUMAN/"10_GPL17586_usable_probe_symbol_mapping.csv")
    by_gene=defaultdict(set)
    by_probe=defaultdict(set)
    for row in probe_rows:
        sym=row["symbol"].upper()
        if sym not in archive_gene:continue
        by_gene[sym].add(row["probe"])
        by_probe[row["probe"]].add(sym)
    must(all(len(s)==1 for s in by_probe.values()),
         "Ambiguous reference gene symbol per source probe")
    tested=[gene for gene in archive_gene if archive_gene[gene]["status"]=="TESTABLE"]
    missing=[gene for gene in archive_gene if archive_gene[gene]["status"]!="TESTABLE"]
    must(len(tested)==27 and len(missing)==2,"Unexpected human panel probe coverage")
    for gene in tested:
        must(len(by_gene[gene])==int(archive_gene[gene]["n_probes"]),
             f"Probe count mapping inconsistent {gene}")
    for gene in missing:
        must(len(by_gene[gene])==0 and archive_gene[gene]["status"]=="NO_UNAMBIGUOUS_MATCH",
             "Unrecognized excluded gene status")
    checksum=read(HUMAN/"10_GEO_input_checksums.csv")
    src=[x for x in checksum if x["file"]=="GSE162955_series_matrix.txt.gz"]
    must(len(src)==1 and src[0]["sha256"]==EXPECTED_SHA256 and
         src[0]["source_url"]==SOURCE,"Source manifest SHA mismatch")
    with tempfile.TemporaryDirectory() as folder:
        path=Path(folder)/"GSE162955_series_matrix.txt.gz"
        request=urllib.request.Request(SOURCE,headers={"User-Agent":"IMACS-Stroke-NVU-Audit/2.0"})
        with urllib.request.urlopen(request,timeout=120) as r, path.open("wb") as h:
            while chunk:=r.read(1024*1024):
                h.write(chunk)
        sha=hashlib.sha256(path.read_bytes()).hexdigest()
        must(sha==EXPECTED_SHA256,f"Published human microarray source SHA changed: {sha}")
        titles,donor_metadata=meta_parse(path)
        raw=source_probe_extract(path,by_probe)
    out_rows=[]
    per_donor=[]
    for gene in archive_gene:
        a=archive_gene[gene]
        if gene in missing:
            out_rows.append(dict(gene=gene,status="NO_UNAMBIGUOUS_PLATFORM_PROBE",
                original_status=a["status"],n_probes=0,n_donors=0,
                source_mean_paired_delta="",source_exact_p="",source_frozen_panel_BH_q="",
                leave_one_out_sign_reversals="",leave_one_out_min_delta="",
                leave_one_out_max_delta="",human_tissue_inference="UNAVAILABLE"))
            continue
        vals={}
        for sample in GSM:
            vals[sample]=statistics.median(raw[probe][sample] for probe in by_gene[gene])
        ds=[]
        for patient,(cl,ic) in PATIENTS.items():
            delta=vals[ic]-vals[cl]
            ds.append(delta)
            saved=original_pairs_by.get((gene,patient))
            must(saved is not None,f"Missing saved human patient row {gene}/{patient}")
            must(close(saved["contralateral_intensity"],vals[cl]) and
                 close(saved["infarct_core_intensity"],vals[ic]) and
                 close(saved["paired_core_minus_contra"],delta),
                 f"Original human probe intensity mismatch {gene}/{patient}")
            per_donor.append(dict(gene=gene,patient_id=patient,
                geo_control_GSM=cl,geo_infarct_GSM=ic,
                paired_core_minus_contra_log2=delta,
                age=donor_metadata.get(ic,{}).get("age",""),
                sex=donor_metadata.get(ic,{}).get("gender",
                       donor_metadata.get(ic,{}).get("sex",""))))
        mean,pval=signflip(ds)
        must(close(mean,a["mean_paired_core_minus_contra_log2_intensity"]) and
             close(statistics.median(ds),a["median_paired_core_minus_contra_log2_intensity"]) and
             close(pval,a["paired_exact_signflip_p_two_sided"]),
             f"Human source-derived statistic differs from archive: {gene}")
        must(int((sum(x>0 for x in ds)))==int(float(a["n_positive_pairs"])) and
             sum(x<0 for x in ds)==int(float(a["n_negative_pairs"])),
             "Patient signed counts inconsistent")
        lom=[statistics.mean(ds[:i]+ds[i+1:]) for i in range(6)]
        reversed_n=sum(x*mean < -1e-12 for x in lom)
        p5=[signflip(ds[:i]+ds[i+1:])[1] for i in range(6)]
        out_rows.append(dict(
            gene=gene,status="TESTABLE",original_status=a["status"],
            n_probes=len(by_gene[gene]),n_donors=6,
            source_mean_paired_delta=mean,source_exact_p=pval,
            source_frozen_panel_BH_q=None,
            leave_one_out_sign_reversals=reversed_n,
            leave_one_out_min_delta=min(lom),
            leave_one_out_max_delta=max(lom),
            leave_one_out_min_exact_p_5pairs=min(p5),
            positive_donors=sum(x>0 for x in ds),
            negative_donors=sum(x<0 for x in ds),
            human_tissue_inference="POSTMORTEM_BULK_CORE_VS_CONTRA_DESCRIPTIVE_NOT_BRAIN_EC"))
    accepted=[r for r in out_rows if r["status"]=="TESTABLE"]
    qvalues=bh([r["source_exact_p"] for r in accepted])
    for row,q in zip(accepted,qvalues):
        must(close(q,archive_gene[row["gene"]]["BH_q_frozen_panel"]),
             f"Reconstructed human BH differs {row['gene']}")
        row["source_frozen_panel_BH_q"]=q
    must(len(per_donor)==162, "Unexpected human source-donor-gene count")
    # Cross-platform directional comparison, not a meta-analysis. Different
    # disease models, baselines, timepoints and compartments. No p-value pooling.
    mouse=read(BASE/"primary"/"06_primary_candidate_exact_permutation_screen.csv")
    endothelial={r["gene"].upper():r for r in mouse if r["celltype"]=="Endothelial"}
    sorted_ec={r["gene"].upper():r for r in read(
        BASE/"independent"/"08_independent_FACS_EC_exact_permutation_4v4.csv")}
    stable=read(BASE/"mouse_author_marker_validation"/
                "16_Zheng_pseudobulk_candidate_annotation_sensitivity.csv")
    stable_ec={r["gene"].upper():r for r in stable if r["celltype"]=="Endothelial"}
    cross=[]
    for r in accepted:
        gene=r["gene"]
        must(gene in endothelial and gene in sorted_ec and gene in stable_ec,
             f"Missing mouse EC contrast {gene}")
        x=float(endothelial[gene]["log2CPM_mean_difference"])
        y=float(sorted_ec[gene]["ischemia_minus_naive_log2CPM"])
        z=float(r["source_mean_paired_delta"])
        stableflag=stable_ec[gene]["sign_concordant"]
        must(stableflag in ("True","False"),"Unexpected EC sensitivity flag")
        cross.append(dict(gene=gene,human_core_vs_contra_mean_log2=z,
            human_exact_six_pair_p=r["source_exact_p"],
            human_BH_q_frozen_panel=r["source_frozen_panel_BH_q"],
            human_leave_one_donor_sign_changes=r["leave_one_out_sign_reversals"],
            mouse_GSE174574_MCAO_vs_sham_24h_log2CPM=x,
            mouse_GSE174574_n_animal_each_arm=3,
            mouse_GSE174574_inhouse_EC_annotation_direction_stable=stableflag,
            mouse_GSE223744_FACS_EC_ischemia_vs_naive_24h_log2CPM=y,
            mouse_GSE223744_n_animal_each_arm=4,
            human_whole_tissue_direction_matches_both_mouse_EC=
                x*y>0 and x*z>0,
            translation_gate="HYPOTHESIS_ONLY_DIFFERENT_SPECIMENS_AND_BASELINES",
            PET_gate="NOT_ESTABLISHED_PROTEIN_ACCESSIBILITY_BINDER_PET_PK"))
    OUT.mkdir(parents=True,exist_ok=True)
    save(OUT/"18_GSE162955_human_source_verified_paired_sensitivity.csv",out_rows)
    save(OUT/"18_GSE162955_human_verified_patient_pair_deltas.csv",per_donor)
    save(OUT/"18_human_mouse_cross_model_direction_descriptive.csv",cross)
    report={
      "human_dataset":"GSE162955", "organism":"Homo sapiens",
      "source_series_sha256":EXPECTED_SHA256,"verified_original_sample_titles":titles,
      "tissue":"postmortem infarct core vs contralateral whole-brain homogenate",
      "pairing":"6 matched donors, exact GSM/title pairs source checked",
      "gene_mapping":"Original notebook GPL17586 unambiguous probes, no new mapping aliases; n27 testable",
      "number_testable_genes":len(accepted), "unmatched_gene_symbols":missing,
      "human_frozen_27gene_FDR_below_0_05":sum(r["source_frozen_panel_BH_q"]<.05 for r in accepted),
      "human_leave_one_patient_sign_changes":{
          r["gene"]:r["leave_one_out_sign_reversals"] for r in accepted
          if r["leave_one_out_sign_reversals"]>0},
      "genes_same_direction_in_human_and_two_mouse_EC_models":[
          r["gene"] for r in cross if r["human_whole_tissue_direction_matches_both_mouse_EC"]],
      "model_comparability":"NO formal cross-species pooled statistics; source controls, tissue compartment, lesion pathology and time differ",
      "human_endothelial_specificity":"NOT MEASURED; homogenate from infarct core (necrosis/cell loss) vs contralateral",
      "statistical_note":"64 exact paired sign-flips; 5-donor leave-one-out has minimum attainable two-sided p=0.0625; panel-only BH (not genome-wide)",
      "protein_PET":"None measured",
    }
    (OUT/"18_human_replication_and_translation_audit.json").write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("HUMAN_REPLICATION_AUDIT PASS",json.dumps({
        "testable":len(accepted),"BH05":report["human_frozen_27gene_FDR_below_0_05"],
        "leave_one_out_sign_flip_genes":len(report["human_leave_one_patient_sign_changes"]),
        "cross_model_match_count":len(report["genes_same_direction_in_human_and_two_mouse_EC_models"]),
        "cross_model_match_genes":report["genes_same_direction_in_human_and_two_mouse_EC_models"],
    }),flush=True)

if __name__=="__main__":
    main()
