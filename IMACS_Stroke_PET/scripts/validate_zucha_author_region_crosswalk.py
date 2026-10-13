#!/usr/bin/env python3
"""Strict original-Zucha Seurat-to-GEO spot-region crosswalk, GSE233814.

No automatic timepoint-based section assignment: there are TWO day-7 sections.
Original region annotations must come from published Mendeley Seurat RDS
DOI:10.17632/gnb2dsjms2.1 (export_zucha_region_annotations.R).
A verified author_sample -> GSM section mapping is required for the four final author sections, or all five GEO source sections when both day-7 sections are present.
No spots are treated as independent animal-level replicates.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import math
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
GEO=ROOT/"results"/"spatial"/"04_spatial_GEO_exact_sample_manifest.csv"
SPOTS=ROOT/"results"/"spatial"/"04_spatial_candidate_values_by_spot_NO_REGIONS.csv.gz"
OUT=ROOT/"results"/"spatial_author_region_qc"
SERIES={
 "GSM7437221": "control",
 "GSM7437222": "day1",
 "GSM7437223": "day3",
 "GSM7437224": "day7",
 "GSM7437225": "day7",
}
EXPECTED_CONDITIONS={
 "control":{"ctrl","control","sham"},
 "day1":{"1dpi","1day","day1","d1"},
 "day3":{"3dpi","3day","day3","d3"},
 "day7":{"7dpi","7day","day7","d7","d7b"},
}

def require(ok,msg):
    if not ok: raise AssertionError(msg)

def read(path):
    with open(path,encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def clean_label(s):
    return "".join(ch for ch in str(s).lower() if ch.isalnum())

def save_csv(p,data):
    require(data,"Refusing empty output")
    with open(p,"w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)

def validate_crosswalk(author,section_map,spot_csv,geo_csv,out,min_spots=100):
    geo=read(geo_csv)
    require(set(r["gsm"] for r in geo)==set(SERIES) and len(geo)==5,
            "Archive must match five original GSE233814 GEO sections")
    require(len(section_map) in (4,5),"Require 4 final author sections or 5 original sections")
    sampled={}
    for r in section_map:
        name,gsm=r["author_sample"],r["gsm"]
        require(bool(name) and gsm in SERIES and name not in sampled,
                "Nonunique or invalid author sample map")
        sampled[name]=gsm
    mapped_gsms=set(sampled.values())
    require(len(mapped_gsms)==len(sampled),"Two author sections cannot map to one GEO GSM")
    require({SERIES[x] for x in mapped_gsms}=={"control","day1","day3","day7"},
            "Original four condition categories need distinct anatomical sections")
    require(len(mapped_gsms.intersection({"GSM7437224","GSM7437225"}))==len(mapped_gsms)-3,
            "Day7 source sections cannot be inferred from timepoint alone")
    require(author and all(("author_barcode" in r and "author_sample" in r and
                             "DetailedRegionAnnoShort" in r and "condition" in r and
                             "original_seurat_spot_id" in r) for r in author),
            "Missing original Seurat sample/spot/region metadata; no inference")
    require(len(set(r["original_seurat_spot_id"] for r in author))==len(author),
            "Duplicated original Seurat spot identifiers")
    aligned={}
    skipped_missing_region=0
    for r in author:
        sample=r["author_sample"]
        require(sample in sampled, "Unmapped author sample ID "+sample)
        gsm=sampled[sample]
        label=clean_label(r["condition"])
        require(label in EXPECTED_CONDITIONS[SERIES[gsm]],
                f"Author timepoint mismatch in manually assigned {sample}->{gsm}")
        bc=r["author_barcode"]
        require(bool(bc),"No original author barcode in "+str(r["original_seurat_spot_id"]))
        key=(gsm,bc)
        require(key not in aligned,"Collision in author GEO/barcode mapping "+str(key))
        region=r["DetailedRegionAnnoShort"].strip()
        if not region:
            skipped_missing_region+=1
            region=None
        aligned[key]=dict(author_region=region,
                          original_condition=r["condition"],
                          original_sample=sample)
    require(len(aligned)>=min_spots*len(sampled),"Insufficient author metadata records")
    seen_spots=set()
    hit_spots=set()
    aggregate=defaultdict(lambda:[0,0.0])
    genes_by_gsm=defaultdict(set)
    with gzip.open(spot_csv,"rt",encoding="utf-8",newline="") as f:
        rdr=csv.DictReader(f)
        require({"gsm","spot_id","gene","log1p_cpm"}.issubset(set(rdr.fieldnames or [])),
                "Unexpected archived spot-expression schema")
        for r in rdr:
            key=(r["gsm"],r["spot_id"])
            require(key[0] in SERIES,"Unrecognized GEO sample in stored spots")
            seen_spots.add(key)
            genes_by_gsm[key[0]].add(r["gene"])
            if key not in aligned or aligned[key]["author_region"] is None:
                continue
            v=float(r["log1p_cpm"])
            require(math.isfinite(v) and v>=0,"Nonfinite archive expression")
            hit_spots.add(key)
            dest=(key[0],aligned[key]["author_region"],r["gene"])
            aggregate[dest][0]+=1
            aggregate[dest][1]+=v
    matched=Counter(gsm for gsm,_ in hit_spots)
    authors=Counter(gsm for gsm,_ in aligned)
    geo_spots=Counter(gsm for gsm,_ in seen_spots)
    report=[]
    for gsm in sorted(mapped_gsms):
        m=matched[gsm]
        require(m>=min_spots,f"Too few exact Seurat/GEO barcode matches {gsm}")
        # The final author Seurat object may be QC-filtered; allow missing
        # GEO spots, but never permit unmatched author spots to be silently used.
        fraction=m/authors[gsm]
        require(fraction>=.90,f"Only {fraction:.1%} of author spots map to GEO: {gsm}")
        report.append(dict(
            gsm=gsm,geo_section=SERIES[gsm],
            author_sample=next(a for a,b in sampled.items() if b==gsm),
            geo_unique_spots=geo_spots[gsm],
            original_author_metadata_spots=authors[gsm],
            exactly_matched_annotated_spots=m,
            fraction_author_spots_matched=fraction,
            distinct_archived_candidate_features=len(genes_by_gsm[gsm]),
            mismatch_gate="PASS_SPOT_IDENTITY_NOT_BIOLOGICAL_INFERENCE"))
    rows=[]
    for (gsm,region,gene),(count,total) in sorted(aggregate.items()):
        rows.append(dict(gsm=gsm,experimental_section=SERIES[gsm],
            author_DetailedRegionAnnoShort=region,gene=gene,
            n_matched_spots=count,mean_spot_log1p_cpm=total/count,
            biological_inference="DESCRIPTIVE_ONLY_ONE_SECTION_NOT_INDEPENDENT_REPLICATES"))
    require(rows,"Original author regions yielded no matched expression")
    out.mkdir(parents=True,exist_ok=True)
    save_csv(out/"17_spatial_original_author_region_match_QC.csv",report)
    save_csv(out/"17_spatial_original_regions_candidate_expression_DESCRIPTIVE.csv",rows)
    summary={
       "dataset":"GSE233814 (subseries of GSE233815)",
       "original_author_RDS":"10.17632/gnb2dsjms2.1",
       "author_sample_to_GSM_mappings_verified_manually":len(sampled),
       "unmapped_GEO_spatial_sections_excluded_from_author_region_transfer":sorted(set(SERIES)-mapped_gsms),
       "original_author_processing_note":"Author 1DP_02 workflow removes D7b after integration because its bregma section differs from others; final author Seurat may legitimately have four, not five, sections.",
       "barcodes_joined_by_exact_GSM_and_author_barcode":True,
       "no_condition_to_GSM_guessing":True,
       "no_spot_level_inferential_pvalues":True,
       "no_core_penumbra_remapping_without_author_field":True,
       "unassigned_original_author_region_rows":skipped_missing_region,
       "sections":report,
       "sample_level_replicates":{
           "control":1,"day1":1,"day3":1,"day7":2},
       "interpretation":"Spatial author region labels descriptively map archived spot RNA, NOT vascular-cell localization or rat/human PET feasibility"
    }
    (out/"17_spatial_region_crosswalk_manifest.json").write_text(
        json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary

def self_test():
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)
        def write(path,rows):
            save_csv(path,rows)
        refs=[dict(gsm=x) for x in SERIES]
        write(p/"geo.csv",refs)
        mapping=[dict(author_sample="section"+str(i),gsm=gsm)
                 for i,gsm in enumerate(SERIES,1)]
        author=[]
        for i,(gsm,condition) in enumerate(SERIES.items(),1):
            for j in range(2):
                author.append(dict(original_seurat_spot_id=f"s{i}_{j}",
                    author_barcode=f"TEST{i}_{j}-1",author_sample=f"section{i}",
                    condition={"control":"Ctrl","day1":"1DPI","day3":"3DPI","day7":"7DPI"}[condition],
                    DetailedRegionAnnoShort=f"author_region_{i}"))
        with gzip.open(p/"spots.gz","wt",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=["gsm","spot_id","gene","log1p_cpm"])
            writer.writeheader()
            for i,gsm in enumerate(SERIES,1):
                for j in range(2):
                    writer.writerow(dict(gsm=gsm,spot_id=f"TEST{i}_{j}-1",
                                         gene="Itm2a",log1p_cpm=str(i)))
        result=validate_crosswalk(author,mapping,p/"spots.gz",p/"geo.csv",
                                  p/"out",min_spots=1)
        require(len(result["sections"])==5,"Missing five-source test sample")
        # The original source paper removed D7b from its final object; a
        # four-section author export must be accepted without inventing a
        # fifth region, provided each retained section maps uniquely.
        result4=validate_crosswalk(author[:8],mapping[:4],p/"spots.gz",
                                   p/"geo.csv",p/"out4",min_spots=1)
        require(len(result4["sections"])==4,"Original final 4-section export rejected")
        require(len(result4["unmapped_GEO_spatial_sections_excluded_from_author_region_transfer"])==1,
                "Missing source section not explicitly excluded")
        # Deliberately collapse both day-7 author_sample ids: must refuse.
        bad=[dict(x) for x in mapping]
        bad[-1]["author_sample"]=bad[-2]["author_sample"]
        try:
            validate_crosswalk(author,bad,p/"spots.gz",p/"geo.csv",
                               p/"bad",min_spots=1)
        except AssertionError:
            pass
        else: raise AssertionError("Could not reject ambiguous day7 mapping")
        print("SPATIAL_REGION_CROSSWALK_SELFTEST PASS: five original or four final author sections; ambiguous day7 rejected")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--author-metadata",type=Path)
    parser.add_argument("--section-map",type=Path)
    a=parser.parse_args()
    if a.self_test:
        self_test()
    else:
        require(a.author_metadata is not None and a.section_map is not None,
                "Original author Seurat export and explicit five-row section map required")
        report=validate_crosswalk(read(a.author_metadata),read(a.section_map),
                                  SPOTS,GEO,OUT)
        print("SPATIAL_AUTHOR_REGION CROSSWALK PASS",json.dumps(report),flush=True)
