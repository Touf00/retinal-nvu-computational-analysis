#!/usr/bin/env python3
"""Source-specific adult human brain FACS-EC GEO manifest discovery.

GSE256490 (Wälchli et al. Nature 2024) is a HUMAN VASCULAR REFERENCE,
not a human acute ischemic stroke dataset. A library is not automatically
a distinct donor. All observed raw GEO labels are preserved for manual QC.
"""
from pathlib import Path
from collections import defaultdict
import csv, gzip, hashlib, json, urllib.request, re, tempfile
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"human_vascular_reference"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256490/soft/GSE256490_family.soft.gz"

def check(ok,msg):
    if not ok: raise AssertionError(msg)

def put(name,rows):
    check(bool(rows),"Empty manifest "+name)
    with (OUT/name).open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"GSE256490_family.soft.gz"
        r=urllib.request.Request(URL,headers={"User-Agent":"IMACS-Stroke-NVU-PET/2.0"})
        with urllib.request.urlopen(r,timeout=180) as h,p.open("wb") as f:
            while item:=h.read(1024*1024): f.write(item)
        sha=hashlib.sha256(p.read_bytes()).hexdigest()
        rows=[];s=None
        with gzip.open(p,"rt",encoding="utf-8",errors="replace") as handle:
            for line in handle:
                line=line.rstrip("\r\n")
                if line.startswith("^SAMPLE = "):
                    if s is not None:rows.append(s)
                    s=dict(gsm=line.split("=",1)[1].strip(),title="",source="",organism="",
                           platform="",characteristics=[],supplementary=[])
                elif line.startswith("^") and s is not None:
                    rows.append(s);s=None
                elif s is not None and line.startswith("!Sample_"):
                    key,_,val=line.partition(" = ")
                    val=val.strip()
                    if key=="!Sample_title":s["title"]=val
                    elif key=="!Sample_source_name_ch1":s["source"]=val
                    elif key=="!Sample_organism_ch1":s["organism"]=val
                    elif key=="!Sample_platform_id":s["platform"]=val
                    elif key.startswith("!Sample_characteristics_ch1"):
                        s["characteristics"].append(val)
                    elif key.startswith("!Sample_supplementary_file"):
                        s["supplementary"].append(val)
        if s is not None:rows.append(s)
    check(50<=len(rows)<=220,"Unexpected GEO sample total "+str(len(rows)))
    check(len(set(s["gsm"] for s in rows))==len(rows),"Duplicate GEO sample IDs")
    check(all(x["organism"]=="Homo sapiens" for x in rows),"Non-human GEO sample")
    allrows=[];candidates=[]
    for s in rows:
        title=s["title"]
        parts=re.split(r"[_\s-]+",title.upper())
        is_tl=("TL" in parts or "TEMPORAL" in parts)
        is_ec=("EC" in parts) and ("UNSORTED" not in parts)
        is_fetal="FETAL" in parts
        # Source paper defines TL as surgically available adult/control
        # temporal-lobe brain material. Neither health nor unique donor IDs
        # is inferred from library name alone.
        proposed=is_tl and is_ec and not is_fetal
        sup=s["supplementary"]
        kinds={}
        for suffix in ("matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz"):
            matches=[x.replace("ftp://","https://",1) for x in sup
                     if x.lower().endswith(suffix)]
            check(len(matches)<=1,f"Duplicate {suffix} {s['gsm']}")
            kinds[suffix]=matches[0] if matches else ""
        record=dict(gsm=s["gsm"],title=title,
            source=s["source"],organism=s["organism"],platform=s["platform"],
            characteristics="; ".join(s["characteristics"]),
            candidate_adult_temporal_lobe_FACS_EC=proposed,
            title_based_donor_id="UNVERIFIED",
            verification_scope="AUTHOR_GEO_METADATA_ONLY",
            matrix_url=kinds["matrix.mtx.gz"],
            features_url=kinds["features.tsv.gz"],
            barcodes_url=kinds["barcodes.tsv.gz"])
        allrows.append(record)
        if proposed:candidates.append(record)
    check(len(candidates)>=2,
          "Too few TL FACS-EC libraries to support donor-level screening")
    missing=[r["gsm"] for r in candidates if not (r["matrix_url"] and
             r["features_url"] and r["barcodes_url"])]
    report={
        "dataset":"GSE256490",
        "paper":"Wälchli et al. Nature 2024, DOI 10.1038/s41586-024-07493-y",
        "source_family_soft":URL,"compressed_source_sha256":sha,
        "total_GEO_libraries":len(rows),
        "candidate_adult_TL_CD31_CD45neg_endothelial_libraries":len(candidates),
        "candidate_libraries_missing_any_matrix_url":missing,
        "candidate_GEO_titles":[{"gsm":r["gsm"],"title":r["title"],
                                "characteristics":r["characteristics"]}
                               for r in candidates],
        "patient_mapping":"DO NOT INFER INDIVIDUAL DONOR COUNTS FROM LIBRARY NAMES",
        "disease_scope":"Non-stroke adult temporal-lobe control vascular reference ONLY",
        "downstream_role":"Healthy cerebral EC target baseline / background expression; not stroke change or PET ligand proof",
        "selection_basis":"GEO title TL and token EC without UNSORTED/FETAL",
    }
    put("19_GSE256490_full_original_GEO_sample_inventory.csv",allrows)
    put("19_GSE256490_adult_TL_FACS_EC_candidate_libraries.csv",candidates)
    (OUT/"19_GSE256490_candidate_inventory_audit.json").write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("HUMAN_VASCULAR_GEO_DISCOVERY PASS",
          json.dumps({"total":len(rows),"candidate_libraries":len(candidates),
                      "missing_source":missing,
                      "candidates":[(r["gsm"],r["title"]) for r in candidates]}),
          flush=True)

if __name__=="__main__":main()
