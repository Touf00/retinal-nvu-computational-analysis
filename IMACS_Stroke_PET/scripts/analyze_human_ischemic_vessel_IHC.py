#!/usr/bin/env python3
"""Pre-registered donor-level analysis for human stroke FFPE vessel IHC.

IMPORTANT: this script contains NO real human staining measurements.
Input is human-generated per-vessel ROI quantitation after blinded QC.
It tests patient-level paired *protein staining*; it CANNOT prove that
any antigen is accessible from the lumen of a living brain microvessel.
Thresholds come from a separate locked assay-control CSV before unblinding.

Usage:
 python IMACS_Stroke_PET/scripts/analyze_human_ischemic_vessel_IHC.py \
   --vessels vessels.csv --thresholds locked_controls.csv --out results/ihc
 python IMACS_Stroke_PET/scripts/analyze_human_ischemic_vessel_IHC.py --self-test
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import random
import statistics
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

TARGETS={"ICAM1","VCAM1"}
REGIONS={"peri_infarct","infarct_core","contralateral","nonstroke_control"}
SEGMENTS={"capillary","postcapillary_venule","arteriole","indeterminate"}
HEADERS={
  "donor_id","block_id","slide_id","vessel_id","target","assay",
  "stain_batch","clinical_group","poststroke_hours","region",
  "segment","endothelial_mask_pixels","target_mean_intensity",
  "local_background_intensity","mask_qc_pass","pathology_region_verified",
  "acquisition_blinded"
}
CONTROL_HEADERS={
 "target","stain_batch","corrected_intensity_cutoff",
 "n_independent_negative_control_fields","positive_control_pass",
 "threshold_locked_before_unblinding"
}
MIN_MICROVESSELS_PER_REGION=10
MIN_DONORS_FOR_EXACT_TEST=6
PILOT_ACUTE_HOURS_MAX=72
RNG_SEED=20261013

def require(ok,msg):
    if not ok: raise ValueError(msg)

def read(path):
    with path.open(newline="",encoding="utf-8") as f:
        rdr=csv.DictReader(f)
        require(rdr.fieldnames is not None,"Empty CSV "+str(path))
        return rdr.fieldnames,list(rdr)

def write(path,rows):
    require(bool(rows),"No output rows for "+str(path))
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def finite(x):
    try:
        value=float(x)
    except (TypeError,ValueError):
        raise ValueError(f"Not numeric: {x}") from None
    require(math.isfinite(value),"Non-finite quantitative intensity")
    return value

def yes(x):
    require(x in {"True","False"},"Expected literal True/False")
    return x=="True"

def exact_signflip(deltas):
    """Paired donor mean difference under 2^n sign assignments, two-sided."""
    n=len(deltas)
    if n<MIN_DONORS_FOR_EXACT_TEST or n>18:
        return None
    observed=statistics.mean(deltas)
    extreme=0
    total=2**n
    for signs in itertools.product((-1,1),repeat=n):
        null=sum(d*s for d,s in zip(deltas,signs))/n
        if abs(null)>=abs(observed)-1e-12:
            extreme+=1
    return extreme/total

def percentile(seq,p):
    s=sorted(seq)
    idx=(len(s)-1)*p
    low=int(math.floor(idx));hi=int(math.ceil(idx))
    return s[low]+(idx-low)*(s[hi]-s[low])

def bootstrap_ci(deltas):
    """Donor-resampling uncertainty interval; descriptive small-N pilot."""
    if len(deltas)<2:
        return None,None
    r=random.Random(RNG_SEED)
    sample=[statistics.mean([deltas[r.randrange(len(deltas))]
                  for _ in range(len(deltas))]) for j in range(3000)]
    return percentile(sample,.025),percentile(sample,.975)

def bh_two(pvals):
    ordered=sorted(pvals,key=lambda k:pvals[k])
    q={}
    prev=1.0
    for rank in range(len(ordered),0,-1):
        k=ordered[rank-1]
        prev=min(prev,pvals[k]*len(ordered)/rank)
        q[k]=prev
    return q

def run(vessel_csv,threshold_csv,out):
    hdr,controls=read(threshold_csv)
    require(CONTROL_HEADERS.issubset(set(hdr)),"Missing locked-control columns")
    thresholds={}
    for c in controls:
        key=(c["target"],c["stain_batch"])
        require(c["target"] in TARGETS and bool(c["stain_batch"]),
                "Unexpected target/batch in assay calibration")
        require(key not in thresholds,"Repeated assay threshold")
        require(int(c["n_independent_negative_control_fields"])>=3,
                "Missing background/isotype control fields")
        require(yes(c["positive_control_pass"]) and
                yes(c["threshold_locked_before_unblinding"]),
                "Assay calibration not independently validated or unlocked")
        val=finite(c["corrected_intensity_cutoff"])
        require(val>0,"Positive threshold must be >0")
        thresholds[key]=val
    require(bool(thresholds),"No independently qualified assay thresholds")

    header,records=read(vessel_csv)
    require(HEADERS.issubset(set(header)),"Missing vessel QC/provenance columns")
    require(records,"No experimental vessel rows")
    seen=set()
    donor_info={}
    accepted=[]
    exclusions=Counter()
    for r in records:
        donor=r["donor_id"]
        require(donor and donor.isascii() and
                not ("@" in donor or " " in donor),
                "Non-anonymized or empty donor ID")
        require(r["target"] in TARGETS,"Unexpected marker")
        require(r["assay"]=="FFPE_IHC",
                "FFPE image table cannot contain live luminal binding data")
        require((r["target"],r["stain_batch"]) in thresholds,
                "No independently locked calibration for marker/batch")
        require(r["clinical_group"] in {"AIS","nonstroke_control"},
                "Unexpected clinical group")
        require(r["region"] in REGIONS and r["segment"] in SEGMENTS,
                "Unrecognized location or vessel phenotype")
        if r["clinical_group"]=="AIS":
            require(r["region"]!="nonstroke_control",
                    "AIS cannot be classified as unrelated donor control")
            hours=finite(r["poststroke_hours"])
            require(hours>0,"Invalid onset to death interval")
        else:
            require(r["region"]=="nonstroke_control" and
                    r["poststroke_hours"]=="",
                    "Control donor assigned ischemic region or time")
            hours=None
        identity=(donor,r["block_id"],r["slide_id"],r["vessel_id"],r["target"])
        require(all(identity),"Missing block/slide/vessel identity")
        require(identity not in seen,"Duplicate per-vessel target/slide ROI")
        seen.add(identity)
        donor_meta=(r["clinical_group"],hours)
        if donor in donor_info:
            require(donor_info[donor]==donor_meta,
                    "Contradictory donor group or stroke interval")
        donor_info[donor]=donor_meta
        area=finite(r["endothelial_mask_pixels"])
        raw=finite(r["target_mean_intensity"])
        bg=finite(r["local_background_intensity"])
        require(area>=0 and raw>=0 and bg>=0,"Negative area/intensity")
        flags=[yes(r[x]) for x in ("mask_qc_pass","pathology_region_verified",
                                  "acquisition_blinded")]
        if not all(flags):
            exclusions["FAILED_MANUAL_BLINDED_QC_OR_UNVERIFIED_REGION"]+=1
            continue
        if area<100:
            exclusions["ENDOTHELIAL_MASK_LT100_PIXELS"]+=1
            continue
        if r["segment"]=="indeterminate":
            exclusions["UNRESOLVED_VESSEL_SEGMENT"]+=1
            continue
        corrected=max(0.0,raw-bg)
        threshold=thresholds[(r["target"],r["stain_batch"])]
        accepted.append(dict(
            donor_id=donor,block_id=r["block_id"],slide_id=r["slide_id"],
            vessel_id=r["vessel_id"],target=r["target"],
            region=r["region"],segment=r["segment"],
            clinical_group=r["clinical_group"],poststroke_hours=hours,
            stain_batch=r["stain_batch"],
            endothelial_mask_pixels=area,
            bg_corrected_endothelial_target_intensity=corrected,
            positive_vessel=int(corrected>threshold),
            prespecified_positive_cutoff=threshold,
            protein_interpretation="FFPE_ENDOTHELIAL_COLOCALIZATION_ONLY_NOT_LIVE_LUMINAL_ACCESS"
        ))
    require(accepted,"All vessels were excluded")
    # Never aggregate until each row has a donor key and a definite vessel identity.
    groups=defaultdict(list)
    for r in accepted:
        groups[(r["donor_id"],r["target"],r["region"],r["segment"])].append(r)
    donor_rows=[]
    for key,values in sorted(groups.items()):
        donor,target,region,segment=key
        donor_rows.append(dict(
            donor_id=donor,target=target,region=region,segment=segment,
            clinical_group=donor_info[donor][0],
            poststroke_hours=donor_info[donor][1],
            n_blinded_QC_pass_vessels=len(values),
            fraction_vessels_protein_positive=
               sum(r["positive_vessel"] for r in values)/len(values),
            mean_corrected_target_intensity=
               statistics.mean(r["bg_corrected_endothelial_target_intensity"]
                               for r in values),
            median_corrected_target_intensity=
               statistics.median(r["bg_corrected_endothelial_target_intensity"]
                                 for r in values),
            evidence_type="PER_DONOR_FFPE_ENDOTHELIAL_SEGMENT_PROTEIN"
        ))
    # Primary analysis is acute (<=72h), CAPILLARY vessel phenotype,
    # infarct border versus contralateral within the same donor.
    # Avoid structural artefacts in necrotic infarct cores and capillary
    # versus postcapillary venule mix. Other contrasts remain descriptive.
    keyed={(r["donor_id"],r["target"],r["region"],r["segment"]):r
           for r in donor_rows}
    paired=[]
    summary=[]
    tests={}
    for target in sorted(TARGETS):
        paired_diffs=[]
        for donor,(group,hours) in sorted(donor_info.items()):
            if group!="AIS" or hours>PILOT_ACUTE_HOURS_MAX:
                continue
            a=keyed.get((donor,target,"peri_infarct","capillary"))
            b=keyed.get((donor,target,"contralateral","capillary"))
            if a is None or b is None:
                continue
            if min(a["n_blinded_QC_pass_vessels"],
                   b["n_blinded_QC_pass_vessels"])<MIN_MICROVESSELS_PER_REGION:
                continue
            delta=(a["fraction_vessels_protein_positive"]-
                   b["fraction_vessels_protein_positive"])
            paired_diffs.append(delta)
            paired.append(dict(donor_id=donor,target=target,
                 acute_hours=hours,segment="capillary",
                 peri_fraction=a["fraction_vessels_protein_positive"],
                 contralateral_fraction=b["fraction_vessels_protein_positive"],
                 paired_fraction_difference=delta,
                 n_peri_vessels=a["n_blinded_QC_pass_vessels"],
                 n_contralateral_vessels=b["n_blinded_QC_pass_vessels"]))
        p=exact_signflip(paired_diffs)
        lo,hi=bootstrap_ci(paired_diffs)
        tests[target]=p
        summary.append(dict(
            target=target,primary_contrast="acute_peri_infarct_minus_contralateral",
            unit_of_inference="PATIENT_DONOR_NOT_VESSEL",
            vessel_segment="capillary",
            n_paired_acute_donors=len(paired_diffs),
            min_paired_microvessels_per_donor_region=MIN_MICROVESSELS_PER_REGION,
            mean_donor_level_positive_fraction_difference=(
                statistics.mean(paired_diffs) if paired_diffs else ""),
            median_donor_level_positive_fraction_difference=(
                statistics.median(paired_diffs) if paired_diffs else ""),
            donor_bootstrap_95CI_low=lo if lo is not None else "",
            donor_bootstrap_95CI_high=hi if hi is not None else "",
            exact_signflip_p=(p if p is not None else ""),
            BH_q_two_targets="",
            status=("DESCRIPTIVE_INSUFFICIENT_DONORS" if p is None else
                    "DONOR_LEVEL_EXACT_PAIRED_PILOT_ONLY"),
            human_luminal_target_accessibility="UNASSESSED_FFPE_CANNOT_TEST",
            PET_target_validation="NOT_TESTED"))
    p_available={x:p for x,p in tests.items() if p is not None}
    q=bh_two(p_available) if len(p_available)==2 else {}
    for row in summary:
        row["BH_q_two_targets"]=q.get(row["target"],"")
    out.mkdir(parents=True,exist_ok=True)
    write(out/"23_qc_passing_ffpe_vessel_rois.csv",accepted)
    write(out/"23_human_ffpe_patient_segment_region_summaries.csv",donor_rows)
    write(out/"23_human_ffpe_per_patient_acute_paired_capillary.csv",paired if paired else [
        {"donor_id":"","target":"","acute_hours":"","segment":"",
        "peri_fraction":"","contralateral_fraction":"",
        "paired_fraction_difference":"","n_peri_vessels":"",
        "n_contralateral_vessels":""}])
    write(out/"23_human_ffpe_donor_level_protein_screen.csv",summary)
    manifest={
      "source_vessels":str(vessel_csv),
      "source_locked_controls":str(threshold_csv),
      "n_experimental_vessel_rows":len(records),
      "n_qc_pass_vessels":len(accepted),
      "excluded_after_predeclared_QC":dict(exclusions),
      "n_distinct_annotated_donors":len(donor_info),
      "primary_endpoint":"Within acute AIS donor peri-infarct vs contralateral CAPILLARY percent IHC-positive vessels; ICAM1 and VCAM1",
      "minimum_acute_paired_patients_for_exact_test":MIN_DONORS_FOR_EXACT_TEST,
      "reproducibility_seed":RNG_SEED,
      "control_thresholds":"Locked in independent batch-specific staining controls before unblinding",
      "inference":"Donor-level paired exact signflip; no vessel-level p-values; two targets BH only if both >=6 pairs",
      "absent_sample_problem":"A missing donor region results in complete-pair exclusion for that target; report missingness, do not impute",
      "limitation":"FFPE protein colocalization does NOT prove surface/luminal accessibility; no PET acquisition",
      "data_privacy":"Only de-identified tissue IDs are accepted, not patient names or source images",
      "biological_result_generated_without_human_data":False
    }
    (out/"23_ihc_preanalysis_audit_manifest.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    return manifest,summary

def selftest():
    """Synthetic test is ephemeral, never archived as clinical results."""
    with tempfile.TemporaryDirectory() as tmp:
        d=Path(tmp)
        c=d/"controls.csv";v=d/"vessels.csv"
        controls=[]
        for target in sorted(TARGETS):
            controls.append(dict(target=target,stain_batch="B1",
                corrected_intensity_cutoff=5,
                n_independent_negative_control_fields=3,
                positive_control_pass="True",
                threshold_locked_before_unblinding="True"))
        write(c,controls)
        rows=[]
        for pid in range(1,7):
            for target in sorted(TARGETS):
                for region in ("peri_infarct","contralateral"):
                    for i in range(12):
                        rows.append(dict(donor_id=f"D{pid:02d}",block_id=region,
                            slide_id=f"S{pid}_{region}",vessel_id=f"V{i}",
                            target=target,assay="FFPE_IHC",stain_batch="B1",
                            clinical_group="AIS",poststroke_hours=24,region=region,
                            segment="capillary",endothelial_mask_pixels=200,
                            target_mean_intensity=12 if region=="peri_infarct" else 4,
                            local_background_intensity=2,mask_qc_pass="True",
                            pathology_region_verified="True",
                            acquisition_blinded="True"))
        write(v,rows)
        report,stats=run(v,c,d/"output")
        require(report["n_distinct_annotated_donors"]==6,"Bad donor count")
        require(all(r["n_paired_acute_donors"]==6 for r in stats),
                "Vessels mistakenly treated as patients")
        require(all(abs(r["exact_signflip_p"]-.03125)<1e-10 for r in stats),
                "Incorrect exact 6-donor signflip")
        require(all(r["human_luminal_target_accessibility"]=="UNASSESSED_FFPE_CANNOT_TEST"
                    for r in stats),"FFPE falsely called luminal access")
        # Missing prespecified batch controls must always hard-fail
        fail_controls=[x for x in controls if x["target"]=="ICAM1"]
        write(c,fail_controls)
        try: run(v,c,d/"invalid")
        except ValueError: pass
        else: raise AssertionError("Failed to refuse uncalibrated marker")
        # Another donor's two regions cannot have contradictory clinical metadata
        write(c,controls)
        wrong=[dict(x) for x in rows]
        wrong[-1]["poststroke_hours"]=48
        write(v,wrong)
        try: run(v,c,d/"invalid2")
        except ValueError:pass
        else:raise AssertionError("Failed inconsistent donor metadata")
        print("HUMAN_FFPE_IHC_PIPELINE_SELFTEST PASS",
              "synthetic test only; no human tissue measurements; n=6 biological patients")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--vessels",type=Path)
    parser.add_argument("--thresholds",type=Path)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--self-test",action="store_true")
    x=parser.parse_args()
    if x.self_test:
        selftest()
    else:
        require(x.vessels and x.thresholds and x.out,
                "--vessels, --thresholds and --out required")
        a,b=run(x.vessels,x.thresholds,x.out)
        print("HUMAN_FFPE_IHC_DONOR_AUDIT PASS",json.dumps({
            "n_donors":a["n_distinct_annotated_donors"],
            "n_qc_vessels":a["n_qc_pass_vessels"],
            "results":b
        }),flush=True)
