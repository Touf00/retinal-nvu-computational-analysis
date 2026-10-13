#!/usr/bin/env python3
"""Fail-closed source/evidence checks for nine IMACS protein/BBB imaging tracks.

Validates internal coherence, provenance and category separation.
This does NOT independently reproduce or confirm literature experiments.
"""
from pathlib import Path
import csv,json
from urllib.parse import urlsplit

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results"/"target_validation_gate"
INPUT=BASE/"22_human_ischemic_vascular_protein_and_tracer_evidence.csv"
OUTPUT=BASE/"22_human_stroke_protein_evidence_consistency.json"
TRACKS={"ICAM1","VCAM1","MCAM","PVR","APLNR","AMOT","SLC3A2","TFRC","MAOB"}

with INPUT.open(encoding="utf-8",newline="") as f:
    rows=list(csv.DictReader(f))
assert len(rows)==9 and {r["gene"] for r in rows}==TRACKS, "Nine target rows missing or duplicated"
assert sum(r["human_stroke_vessel_protein"]=="DIRECT_HUMAN_AUTOPSY" for r in rows)==2, (
    "Human-autopsy status changed; review original source before altering gate"
)
for row in rows:
    gene=row["gene"]
    if gene in {"ICAM1","VCAM1"}:
        assert row["human_stroke_vessel_protein"]=="DIRECT_HUMAN_AUTOPSY"
    elif gene=="MAOB":
        assert row["human_stroke_vessel_protein"]=="NOT_VESSEL_SURFACE_TARGET"
    else:
        assert row["human_stroke_vessel_protein"]=="NOT_IDENTIFIED_IN_REVIEWED_HUMAN_STROKE_BRAIN_TISSUE"
    assert "GO" not in row["decision"] or "NOT_STROKE_PET_GO" in row["decision"], (
        "No human ischemic brain-vessel target has achieved a full PET gate"
    )
    assert "PET" in row["in_vivo_PET_validation"] or gene in {"ICAM1","VCAM1","MCAM","PVR"}
    for field in ("primary_url","secondary_url"):
        u=urlsplit(row[field])
        assert u.scheme=="https" and bool(u.netloc), f"Missing source {gene}: {field}"
    if row["qualifying_url"]:
        u=urlsplit(row["qualifying_url"])
        assert u.scheme=="https" and bool(u.netloc), f"Invalid third source {gene}"
    if gene=="MAOB":
        assert "INTRACELLULAR" in row["intravascular_epitope_status"]
    if gene=="AMOT":
        assert "sCD146" in row["in_vivo_PET_validation"]
    if gene=="MCAM":
        assert "subcutaneous" in row["in_vivo_PET_validation"]

report={
 "n_candidates":len(rows),
 "genes":sorted(TRACKS),
 "direct_human_ischemic_vascular_IHC":[r["gene"] for r in rows if r["human_stroke_vessel_protein"]=="DIRECT_HUMAN_AUTOPSY"],
 "direct_human_ischemic_vascular_PET_verified":[],
 "human_brain_BBB_protein_nonstroke":[r["gene"] for r in rows if r["gene"] in {"TFRC","SLC3A2"}],
 "preclinical_ischemic_cerebral_binding_or_protein":["ICAM1","VCAM1","MCAM","PVR"],
 "protein_surface_access_human_ischemia":"NOT_PROVEN_FOR_ANY_TARGET",
 "documented_stroke_parenchymal_PET":["MAOB-context THK5351 (crossbinding; not specific vascular MAOB)"],
 "literature_review_is_not_original_protein_experiment":True,
 "purpose":"Evidence gate/citation consistency only, not PET clearance or full systematic literature review"
}
OUTPUT.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
print("HUMAN_STROKE_PROTEIN_EVIDENCE_GATE PASS",json.dumps(report),flush=True)
