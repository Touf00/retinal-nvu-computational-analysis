#!/usr/bin/env python3
"""Make a 29-gene cross-species TRANSLATION GATE, not a PET ranking.

Uses only previously executed, source-controlled and archived datasets:
GSE174574 mouse EC 3+3 sham, GSE223744 sorted mouse EC 4+4 naive,
GSE225948 day2 pooled EC 4+4 (batch confounded),
GSE162955 human postmortem paired 6 (bulk infarct core),
GSE256490 human adult surgical-control FACS EC 6 (no stroke).
Never pools effect estimates or reuses cells as biological replicates.
"""
from pathlib import Path
import csv,json,math
from collections import Counter,defaultdict

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results"
OUT=BASE/"target_validation_gate"
OUT.mkdir(parents=True,exist_ok=True)

def read(path):
    with path.open(encoding="utf-8",newline="") as h:
        return list(csv.DictReader(h))

def key(rows,column):
    d={}
    for r in rows:
        k=r[column].upper()
        if k in d:raise AssertionError("Duplicated gene "+k)
        d[k]=r
    return d

def need(ok,msg):
    if not ok:raise AssertionError(msg)

def render(name,items):
    need(len(items)>0,"No gate output")
    with (OUT/name).open("w",encoding="utf-8",newline="") as h:
        wr=csv.DictWriter(h,fieldnames=list(items[0]),extrasaction="ignore")
        wr.writeheader();wr.writerows(items)

def numeric(v):
    if v is None or str(v)=="":
        return None
    out=float(v)
    need(math.isfinite(out),"Nonfinite source effect")
    return out

def main():
    original=read(ROOT/"config"/"candidate_catalog.csv")
    notes=defaultdict(list)
    for r in original:
        if r["role"] not in notes[r["gene"].upper()]:
            notes[r["gene"].upper()].append(r["role"])
    genes=list(dict.fromkeys(r["gene"].upper() for r in original))
    need(len(genes)==29,"Expected frozen 29 genes, no addition")
    acute=key([x for x in read(BASE/"primary"/"06_primary_candidate_exact_permutation_screen.csv")
               if x["celltype"]=="Endothelial"],"gene")
    independent=key(read(BASE/"independent"/"08_independent_FACS_EC_exact_permutation_4v4.csv"),"gene")
    stable=key([r for r in read(BASE/"mouse_author_marker_validation"/"16_Zheng_pseudobulk_candidate_annotation_sensitivity.csv")
                if r["celltype"]=="Endothelial"],"gene")
    pooled=key(read(BASE/"independent_day2"/"11_GSE225948_day2_EC_exact_permutation.csv"),"gene")
    human=key(read(BASE/"human_replication_qc"/"18_GSE162955_human_source_verified_paired_sensitivity.csv"),"gene")
    human_ec=key(read(BASE/"human_vascular_reference"/"20_GSE256490_human_TL_candidate_detection_across_six_preparations.csv"),"gene")
    need(set(acute)==set(independent)==set(stable)==set(human)==set(human_ec)==set(genes),
         "Source did not cover full frozen candidate panel")
    need(len(pooled)==28,"Pooled EC frozen target coverage changed")
    mainrows=[]
    for g in genes:
        a=acute[g];s=independent[g];v=stable[g];h=human[g];e=human_ec[g];p=pooled.get(g)
        x=numeric(a["log2CPM_mean_difference"])
        y=numeric(s["ischemia_minus_naive_log2CPM"])
        z=numeric(h["source_mean_paired_delta"])
        baseline=numeric(e["median_fraction_cells_detected_across_covered_sources"])
        need(v["sign_concordant"]=="True","Endothelial original mouse reannotation direction changed: "+g)
        mouseagree=(x*y>0)
        humanagree=(mouseagree and z is not None and x*z>0)
        qhuman=numeric(h["source_frozen_panel_BH_q"])
        qmouse=numeric(s["BH_q_across_frozen_gene_panel"])
        need(qhuman is None or qhuman>=.05,
             "Unexpected significant human result: inspect separately")
        pooled_effect=numeric(p["stroke_minus_sham_mean_log1p_10000"]) if p else None
        row=dict(
            gene=g,
            frozen_catalog_roles="; ".join(notes[g]),
            human_adult_surgical_control_TL_EC_sources=6,
            human_adult_EC_symbol_covered_preparations=int(e["n_source_symbol_covered"]),
            human_adult_EC_median_fraction_transcript_detected=baseline,
            human_adult_EC_n_preps_with_ge_10pct_detection=int(
                e["n_preparations_with_detection_ge_10pct"]),
            human_adult_EC_measurement="FACS_SORTED_CD31+CD45-_SURGICAL_CONTROL_NOT_STROKE",
            mouse_acute_MCAO_vs_sham_log2CPM=x,
            mouse_acute_n_independent_libraries_each_arm=3,
            mouse_acute_exact_p_min_possible=0.1,
            mouse_acute_EC_reannotation_same_direction=True,
            mouse_sorted_FACS_EC_ischemia_vs_naive_log2CPM=y,
            mouse_sorted_EC_exact_p=numeric(s["permutation_two_sided_p"]),
            mouse_sorted_EC_panel_BH_q=qmouse,
            mouse_sorted_EC_naive_not_sham=True,
            both_mouse_EC_contrasts_same_direction=mouseagree,
            mouse_GSE225948_pooled_4v4_day2_EC_log_mean_effect=pooled_effect,
            mouse_GSE225948_pooled_preparation_batch_confounded=True,
            human_six_donors_core_vs_contralateral_mean_log2=z,
            human_whole_tissue_paired_exact_p=numeric(h["source_exact_p"]),
            human_whole_tissue_panel_BH_q=qhuman,
            human_whole_tissue_source_probe_testable=(h["status"]=="TESTABLE"),
            human_whole_tissue_loo_sign_reversals=(
                int(h["leave_one_out_sign_reversals"]) if h["status"]=="TESTABLE" else ""),
            human_whole_tissue_matches_both_mouse_EC_directions=humanagree,
            convergent_direction_is_not_validation=True,
            model_cross_species_stats="NOT_POOLED_NOT_META_ANALYZED",
            human_celltype_specific_stroke_effect="NOT_AVAILABLE",
            protein_surface_or_abluminal_access="UNVERIFIED",
            PET_binder_affinity_specificity="UNVERIFIED",
            normal_organ_background_and_BBB_leak_confounds="UNVERIFIED",
            PET_TARGET_GO="NO",
        )
        mainrows.append(row)
    render("21_frozen_29_gene_human_mouse_evidence_gates.csv",mainrows)
    by={x["gene"]:x for x in mainrows}
    # Fixed hypothesis-driven *measurement* shortlist. It does not claim that
    # every named protein is luminal or PET-targetable.
    hypotheses=[
      ("MCAM","surface-vessel/mural localization hypothesis","CD146 different from sCD146 ligand AMOT; human EC baseline variable","capillary/venule vs pericyte/VSMC protein, normal-organs, ligand cross-reactivity"),
      ("ICAM1","vascular inflammation surface hypothesis","broad organ endothelial activation, human core statistics not significant","human brain ischemia EC surface IHC and specific blood-pool/BBB-corrected binder"),
      ("PVR","immune-endothelial interface hypothesis","human baseline EC low to intermediate; no validated brain PET binder here","human endothelial surface protein and nonvascular/background binding"),
      ("VCAM1","vascular inflammation comparator hypothesis","human postmortem bulk direction differs from mouse acute; normal systemic inflammation","species-specific brain EC surface protein, timepoint and organ sink"),
      ("APLNR","host APJ falsification hypothesis","mouse acute and sorted EC direction negative, human normal EC RNA sparse","luminal APJ protein after ischemia, receptor blocking vs ischemic blood delivery"),
      ("AMOT","host sCD146 binding-target falsification hypothesis","human normal EC RNA low; AMOT surface fraction and accessibility unresolved","human brain ischemic endothelial accessible AMOT epitope, sCD146 competition and autoradiography"),
      ("SLC3A2","BBB delivery platform comparator","strong human normal EC RNA; transporter expression is NOT intact cargo transcytosis","construct-specific exposure, capillary depletion and target engagement"),
      ("TFRC","BBB delivery platform comparator","substantial normal human EC expression but mouse ischemic RNA decreases","TfR shuttle construct binding, capillary depletion and drug effect"),
      ("MAOB","astrocytic PET biology comparator","human EC expression limited; astrocyte modality, not vascular-surface marker","reactive astrocyte-specific enzyme/protein activity and tracer specificity")
    ]
    lines=[
      "# IMACS focused protein-validation hypotheses — NOT PET nominations",
      "",
      "Frozen panel: all 29 genes remain in the machine-readable evidence table. The following nine are hypothesis-focused laboratory comparators (six vascular imaging hypotheses, two delivery handles, one astrocyte target), not genome-wide statistically prioritized winners.",
      "",
      "Source: six GSE174574 mouse 24h MCAO/Sham libraries; eight GSE223744 FACS sorted mouse ischemia/naive libraries; GSE225948 pooled stroke/sham 4+4 with batch concerns; six GSE162955 paired HUMAN postmortem whole-tissue donors; SIX independently labeled human adult temporal-lobe control FACS brain EC samples from GSE256490. Different estimands MUST NOT be meta-analyzed.",
      "",
      "**Human six-paired infarct-core whole-tissue dataset has ZERO genes with q<0.05. Mouse exact 3v3 permutation p-values have floor 0.10. No candidate is PET-ready.**",
      "",
      "| Gene | Track | Source-aware caution | Next discriminating experiment |",
      "|---|---|---|---|"
    ]
    for gene,track,caution,experiment in hypotheses:
        lines.append(f"| {gene} | {track} | {caution} | {experiment} |")
    lines.extend([
       "",
       "## Evidence flags for the nominated protein-hypothesis set",
       "",
       "| Gene | Human normal EC median RNA detection | Mouse EC 2-cohort direction | Human whole tissue 6-pair mean direction | Human donor leave-one-out sign changes |",
       "|---|---:|---|---|---:|"
    ])
    for gene,*_ in hypotheses:
        r=by[gene]
        md=r["human_adult_EC_median_fraction_transcript_detected"]
        joint="Concordant" if r["both_mouse_EC_contrasts_same_direction"] else "Discordant"
        core=r["human_six_donors_core_vs_contralateral_mean_log2"]
        core_s=("increase" if core>0 else "decrease") if core is not None else "not mappable"
        lines.append(f"| {gene} | {md:.2%} | {joint} | {core_s} | {r['human_whole_tissue_loo_sign_reversals']} |")
    lines.extend([
      "",
      "## Scientific decision gate",
      "",
      "- Low transcript detection in control human EC may reflect true low baseline, sampling/capture dropout, anatomical environment or regulation. It cannot prove protein absent or rule out induction after stroke.",
      "- FACS-sorted human controls are surgical temporal-lobe samples, including epilepsy-referenced filenames; not baseline for all brain regions nor stroke controls.",
      "- Human stroke six-donor core vs contralateral is whole-brain postmortem tissue, affected by cell death and heterogeneous composition. No endothelial-specific paired ischemic effect is established.",
      "- GSE225948 EC2 fraction is highly shifted in pooled stroke vs sham preparations with unexchangeable batch dates; not independent disease attribution.",
      "- Protein localization, accessible epitope, blocking, competition, radiochemistry, physiological transport, normal-organ sink, BBB leak and perfusion confounding are mandatory separate stages.",
      "- **Only after a small protein-screen shortlist survives source- and vessel-validated ischemic localization should the original Zucha spatial DetailedRegionAnnoShort metadata be used to test region specificity.**",
      "",
      "Method: no statistical meta-analysis, no p-value combination, no whole-tissue-to-cell-type causal claims. Human tissue direction matching mouse EC is exploratory only."
    ])
    (OUT/"21_focused_protein_validation_hypotheses.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    report=dict(panel_n=len(mainrows),
      mouse_two_EC_cohorts_direction_concordant=sum(r["both_mouse_EC_contrasts_same_direction"] for r in mainrows),
      mapped_human_paired_genes=sum(r["human_whole_tissue_source_probe_testable"] for r in mainrows),
      mouse_two_EC_cohorts_plus_human_bulk_direction_concordant=sum(r["human_whole_tissue_matches_both_mouse_EC_directions"] for r in mainrows),
      human_core_BH_significant=0,
      experimentally_focus_genes=[g for g,*_ in hypotheses],
      experiment_shortlist_only=True,
      PET_target_GO_n=0,
      note="Human healthy adult surgical control endothelial expression does not validate stroke upregulation or PET binding")
    (OUT/"21_adjudication_summary.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("INTEGRATED_TRANSLATION_GATE PASS",json.dumps(report),flush=True)

if __name__=="__main__":main()
