#!/usr/bin/env python3
"""Audit result provenance and statistical claims in the IMACS analysis branch.

This checks analysis outputs, not PET target validity. It is run automatically
after successful workflow artifacts have been transferred to the branch.
"""
from pathlib import Path
import json
import pandas as pd

BASE = Path("IMACS_Stroke_PET/results")
checks = []
def add(name, status, detail):
    checks.append({"name":name, "status":status, "detail":str(detail)})

def csv(path):
    return pd.read_csv(BASE / path)

try:
    qc=csv("primary/02_QC_audit.csv")
    add("mouse_acute_samples","PASS" if qc.sample_id.nunique()==6 else "FAIL",f"{len(qc)} sample rows")
    x=csv("primary/06_primary_candidate_exact_permutation_screen.csv")
    add("mouse_acute_animal_level_permutations",
        "PASS" if x.n_animals_per_group.eq(3).all() and x.exact_label_permutation_two_sided_p.min()>=.099999 else "FAIL",
        f"{len(x)} candidate x compartment rows; minimum exact p={x.exact_label_permutation_two_sided_p.min():.3f}; cannot claim p<0.05")
    failures=csv("primary/06_DESeq2_failure_audit.csv")
    full_de_available = all((BASE/f"primary/06_DESeq2_{ct}.csv").exists()
                            for ct in ["Endothelial","Astrocyte","Pericyte","VSMC"])
    de_state = ("PASS" if len(failures)==0 and full_de_available else "NOT_COMPLETE")
    add("mouse_acute_DESeq2",de_state,
        f"{len(failures)} unresolved errors; complete four-compartment analysis: {full_de_available}; "
        "n=3/group parametric Wald p-values require careful scrutiny, exact permutation floor 0.10")
except (FileNotFoundError, KeyError) as e:
    add("mouse_acute","NOT_COMPLETE",str(e))

try:
    t=csv("temporal/02_temporal_cell_counts.csv")
    add("mouse_temporal","PASS" if len(t)==4 else "CHECK",f"{len(t)} deposited conditions; ONE SAMPLE per condition, descriptive only")
    cols=csv("temporal/02_candidate_timecourse_descriptive.csv").columns
    add("mouse_temporal_pseudoreplication","PASS" if not {"pvalue","padj","fdr"}.intersection(cols) else "FAIL",
        "No spurious cell-level inferential p-values")
except (FileNotFoundError,KeyError) as e:
    add("mouse_temporal","NOT_COMPLETE",str(e))

try:
    design=csv("rat/03_rat_GEO_design_manifest.csv")
    qc=csv("rat/03_rat_48h_QC_audit.csv")
    add("rat_48h","PASS" if len(design)==6 and len(qc)==6 else "FAIL",
        f"{len(design)} condition/hemisphere matrices; individual rat identities NOT recoverable")
    expr=csv("rat/03_rat_48h_candidate_expression_by_group_and_hemisphere.csv")
    low=expr.loc[expr.celltype.eq("Endothelial"),["prefix","n_cells"]].drop_duplicates()
    add("rat_endothelium_QC","CHECK" if low.n_cells.min()<50 else "PASS",
        f"Minimum endothelial-labelled nuclei in a group = {low.n_cells.min()}, NOT suitable for vascular validation without improved annotation")
    missing=expr.loc[expr.mean_log1p.isna(),"gene"].nunique()
    add("rat_missing_feature_symbols","CHECK" if missing else "PASS",
        f"{missing} of frozen genes could not be matched to reference feature symbols; unmatchable is NOT biological absence")
except (FileNotFoundError,KeyError) as e:
    add("rat_48h","NOT_COMPLETE",str(e))

try:
    sam=csv("spatial/04_spatial_candidate_expression_by_section.csv")
    add("spatial_sections","PASS" if sam.gsm.nunique()==5 else "FAIL",
        f"{sam.gsm.nunique()} Visium sections: Ctrl x1, day1 x1, day3 x1, day7 x2")
    add("spatial_celltypes","CHECK","Mixed Visium spots; no core/periphery annotations; do not assign cell type or lesion location")
except (FileNotFoundError,KeyError) as e:
    add("spatial_sections","NOT_COMPLETE",str(e))

try:
    facs=csv("independent/08_independent_FACS_EC_exact_permutation_4v4.csv")
    test=facs[facs.status=="TESTABLE"]
    add("independent_sorted_endothelial_4v4",
        "PASS" if len(test)>0 and test.n_ischemia.eq(4).all() and test.n_naive.eq(4).all()
          and test.permutation_two_sided_p.min()>=2/70-1e-12 else "FAIL",
        f"{len(test)} frozen genes, 4 ischemia vs 4 naive control; permutation p min={test.permutation_two_sided_p.min():.4f}; naïve controls not shams")
    astro=csv("independent/08_independent_astro_endfoot_paired_3animals.csv")
    val=astro[astro.status=="TESTABLE"]
    add("independent_astro_endfeet_3paired",
        "PASS" if len(val)>0 and val.n_pairs.eq(3).all() else "CHECK",
        f"{len(val)} genes with 3 contralateral/ipsilateral pairs; no population-level significance")
except (FileNotFoundError,KeyError) as e:
    add("independent_vascular_astro","NOT_COMPLETE",str(e))

try:
    additional=csv("spatial_day3/09_day3_spatial_frozen_panel_descriptive.csv")
    m=(BASE/"spatial_day3/09_spatial_run_manifest.json")
    meta=json.loads(m.read_text())
    add("new_2026_independent_spatial_day3",
        "PASS" if meta.get("n_barcodes",0)>1000 and not meta.get("target_protein_or_PET_validation",False) else "CHECK",
        f"GSE337578 one library; {meta.get('n_barcodes')} spatial units; no sham control or validated region masks; descriptive only")
except (FileNotFoundError,KeyError,ValueError) as e:
    add("new_2026_independent_spatial_day3","NOT_COMPLETE",str(e))

try:
    human=csv("human/10_human_paired_core_vs_contralateral_candidate_screen.csv")
    tested=human.loc[human.status.eq("TESTABLE")]
    valid=not tested.empty and tested.n_pairs.eq(6).all() and tested.paired_exact_signflip_p_two_sided.min()>=2/64-1e-12
    add("human_six_paired_brain_infarct",
        "PASS" if valid else "CHECK",
        f"6 matched postmortem brain tissue donors, {len(tested)} mapped candidate genes. "
        "Microarray core-vs-contralateral; not vessel-specific and not PET proof.")
except (FileNotFoundError,KeyError,ValueError) as e:
    add("human_six_paired_brain_infarct","NOT_COMPLETE",str(e))

try:
    day2=csv("independent_day2/11_GSE225948_day2_EC_exact_permutation.csv")
    subtype=csv("independent_day2/11_within_author_EC_subtype_sensitivity.csv")
    ec1=subtype.loc[subtype.author_subtype.eq("EC1") & subtype.status.eq("TESTABLE_4v4_pooled_preps")]
    ec2=subtype.loc[subtype.author_subtype.eq("EC2")]
    safe=(not day2.empty and day2.n_sham_preparations.eq(4).all()
          and day2.n_stroke_preparations.eq(4).all()
          and day2.exact_two_sided_p.min()>=2/70-1e-12
          and len(ec1)==len(day2) and len(ec2)==len(day2)
          and ec2.status.eq("NO_INFERENCE_low_or_missing_group_sample").all())
    qsig=int((day2.BH_q_across_frozen_panel<.05).sum())
    subq=int((ec1.BH_q_within_subtype_frozen_panel<.05).sum())
    add("GSE225948_day2_sham_vs_stroke_pooled_FACS_EC",
        "PASS" if safe else "FAIL",
        f"{len(day2)} testable frozen genes; {qsig} whole-EC and {subq} within-EC1 genes at BH q<.05; "
        "n4+4 independent pooled preparations, NOT individual animals; "
        "non-paired processing dates/batches, author EC-subtype composition shift; source counts fractional")
except (FileNotFoundError,KeyError,ValueError) as e:
    add("GSE225948_day2_FACS_EC","NOT_COMPLETE",str(e))

try:
    late=csv("independent_late/12_GSE225948_day2_EC_exact_permutation.csv")
    sublate=csv("independent_late/12_within_author_EC_subtype_sensitivity.csv")
    late1=sublate.loc[sublate.author_subtype.eq("EC1") & sublate.status.eq("TESTABLE_4v4_pooled_preps")]
    safe=(not late.empty and late.n_day2_stroke_preparations.eq(4).all()
        and late.n_day14_stroke_preparations.eq(4).all()
        and late.exact_two_sided_p.min()>=2/70-1e-12 and len(late1)==len(late))
    add("GSE225948_stroke_day14_vs_day2_pooled_FACS_EC",
        "PASS" if safe else "FAIL",
        f"{len(late)} genes, {(late.BH_q_across_frozen_panel<.05).sum()} whole-EC and "
        f"{(late1.BH_q_within_subtype_frozen_panel<.05).sum()} within-EC1 BH discoveries; "
        "4 pooled stroke preparations/timepoint, UNPAIRED, no sham day14 control; "
        "age/sex matched by source labels, processing batch and subtype composition unresolved")
except (FileNotFoundError,KeyError,ValueError) as e:
    add("GSE225948_stroke_day14_vs_day2","NOT_COMPLETE",str(e))

BASE.mkdir(exist_ok=True,parents=True)
(BASE/"scientific_audit.json").write_text(json.dumps(checks,indent=2))
lines=["# Automated /TEAMRED analysis audit","",
       "Generated from committed processed results. Execution success is not scientific validation.","",
       "| Check | Status | Evidence / limitation |","|---|---|---|"]
for a in checks:
    lines.append("| "+a["name"]+" | "+a["status"]+" | "+a["detail"].replace("|","/").replace("\n"," ")+" |")
lines += ["", "## Mandatory interpretation constraints",
    "- Do not pool cells/spots as biological replicates.",
    "- Mouse GSE174574: n=3/group, minimum attainable exact two-sided permutation p=0.1.",
    "- Temporal GSE227651: one sample/condition, trajectories not statistically replicated.",
    "- Rat GSE250245: six group/hemisphere matrices, not independent individual-animal transcriptomes.",
    "- Spatial GSE233814: mixed capture spots and unavailable lesion-region labels in raw count matrices.",
    "- Rat annotation and missing feature symbols require additional review before vascular validation.",
    "- GSE225948: source matrices are fractional; n4+4 biological pooled preparations, not individual mice.",
    "- Sham and stroke day2 samples are unpaired batches; EC subtype composition differs sharply.",
    "- Day14 versus day2 compares different stroke-only preparations; no sham at day14.",
    "- RNA expression cannot establish protein, PET binding, BBB transfer, or human translation.",
    "", "## Recovery", "",
    "All outputs are inside `IMACS_Stroke_PET/results/<dataset>/` on the isolated branch.",
    "Each dataset's `source_run_id.txt` identifies the successful GitHub Actions execution.",
    "Raw sequencing inputs are downloaded from GEO and have checksums in per-dataset outputs.",
]
(BASE/"SCIENTIFIC_AUDIT.md").write_text("\n".join(lines)+"\n")
for a in checks:print(a["name"],a["status"],a["detail"])
