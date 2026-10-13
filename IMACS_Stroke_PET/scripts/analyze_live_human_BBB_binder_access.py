#!/usr/bin/env python3
"""FFPE-independent HUMAN BBB-on-chip live luminal binding feasibility audit.

This is NOT a human stroke assay or a PET analysis. Measurements must come
from living nonpermeabilized human BBB lumen preparations. A biological
replicate is an independent cell donor/differentiation batch; chips and images
are technical units nested within a donor. The full workflow is a prospective
template, smoke-tested only with synthetic values.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics,tempfile
from collections import defaultdict
from pathlib import Path

ARMS={"target_ligand","size_matched_nonbinder","competitor_blocked_target"}
TARGETS={"ICAM1","VCAM1"}
CONDITIONS={"baseline","TNF_activation","OGD_reoxygenation"}
FIELDS={"donor_batch_id","technical_set_id","chip_id","target","condition",
"assay_arm","assay","label_batch","ligand_size_kDa","input_probe_AUC",
"luminal_retained_probe_signal","luminal_endothelial_area","extravasated_probe_signal",
"permeability_70kDa","viability_pass","barrier_structure_pass",
"probe_added_from_lumen","nonpermeabilized","images_blinded"}
CTRL={"label_batch","max_acceptable_dextran_70kDa_permeability",
"barrier_limit_locked_before_condition_unblinding","positive_live_surface_control_pass"}

def require(value,message):
    if not value:raise ValueError(message)

def number(s):
    n=float(s)
    require(math.isfinite(n),"Nonfinite instrument measurement")
    return n

def boolean(s):
    require(s in ("True","False"),"Invalid boolean flag")
    return s=="True"

def read(path):
    with path.open(newline="",encoding="utf-8") as h:
        x=csv.DictReader(h)
        require(x.fieldnames is not None,"No CSV header")
        return set(x.fieldnames),list(x)

def write(path,values):
    require(values,"No rows to write")
    with path.open("w",newline="",encoding="utf-8") as h:
        wr=csv.DictWriter(h,fieldnames=list(values[0]),extrasaction="ignore")
        wr.writeheader();wr.writerows(values)

def run(chip_csv,locked_controls_csv,destination):
    cfields,controls=read(locked_controls_csv)
    require(CTRL.issubset(cfields),"Missing locked barrier controls")
    control={}
    for c in controls:
        batch=c["label_batch"]
        require(batch and batch not in control,"Duplicate/empty calibration batch")
        require(boolean(c["barrier_limit_locked_before_condition_unblinding"]) and
                boolean(c["positive_live_surface_control_pass"]),
                "Barrier or binding control was not independently qualified")
        threshold=number(c["max_acceptable_dextran_70kDa_permeability"])
        require(threshold>0,"Non-positive permeability upper limit")
        control[batch]=threshold
    hdr,rows=read(chip_csv)
    require(FIELDS.issubset(hdr),"Missing live lumen controls or chip provenance")
    require(rows,"No live endothelial model data")
    sets=defaultdict(dict)
    global_chips=set()
    for r in rows:
        require(r["target"] in TARGETS and r["condition"] in CONDITIONS,
                "Unknown antigen or experiment")
        require(r["assay"]=="LIVE_LUMINAL_FLOW","FFPE may not substitute for live luminal assay")
        require(all(boolean(r[c]) for c in (
                "viability_pass","barrier_structure_pass",
                "probe_added_from_lumen","nonpermeabilized","images_blinded")),
                "Failed essential viability, luminal or blind imaging QC")
        require(r["donor_batch_id"] and r["technical_set_id"] and r["chip_id"],
                "Missing biological or technical unit")
        require(r["assay_arm"] in ARMS,"Missing nonbinding or competition arm")
        require(r["label_batch"] in control,"Unknown or unlocked permeability calibration")
        require(r["chip_id"] not in global_chips,"One chip reused across independent assay arms")
        global_chips.add(r["chip_id"])
        p=number(r["permeability_70kDa"])
        require(p>=0,"Invalid negative permeability")
        # Fail rather than interpreting leakage as target-specific binding.
        require(p<=control[r["label_batch"]],
                "Barrier permeability exceeded locked acceptable maximum")
        size=number(r["ligand_size_kDa"])
        auc=number(r["input_probe_AUC"])
        area=number(r["luminal_endothelial_area"])
        ret=number(r["luminal_retained_probe_signal"])
        extra=number(r["extravasated_probe_signal"])
        require(size>0 and auc>0 and area>0 and ret>=0 and extra>=0,
                "Negative or zero probe/boundary measurement")
        key=(r["donor_batch_id"],r["target"],r["condition"],r["technical_set_id"])
        require(r["assay_arm"] not in sets[key],"Duplicate assay arm within set")
        sets[key][r["assay_arm"]]=dict(
            auc=auc,area=area,signal=ret,extra=extra,
            size=size,permeability=p,label_batch=r["label_batch"])
    experimental=[]
    for (donor,target,condition,set_id),arms in sorted(sets.items()):
        require(set(arms)==ARMS,f"Each set requires target/nonspecific/blocking in {set_id}")
        a,b,c=(arms[x] for x in ("target_ligand","size_matched_nonbinder",
                                 "competitor_blocked_target"))
        require(len({x["label_batch"] for x in arms.values()})==1,
                "Calibration/stain batches differ between matched ligands")
        require(max(x["size"] for x in arms.values())/
                min(x["size"] for x in arms.values())<=1.2,
                "Nonbinding ligand not molecular-size matched")
        require(max(x["auc"] for x in arms.values())/
                min(x["auc"] for x in arms.values())<=1.2,
                "Dose/exposure must be matched; raw SUV cannot test binder specificity")
        # Reporter normalizes exposure AND endothelial luminal area.
        def signal(x):return x["signal"]/x["area"]/x["auc"]
        targeted=signal(a)
        nonspecific=signal(b)
        blocked=signal(c)
        net=targeted-nonspecific
        # Values are descriptive. Off-target baselines and target engagement
        # require antigen knockout/competition by additional orthogonal tests.
        competitor_reduction=(targeted-blocked)/targeted if targeted>0 else None
        experimental.append(dict(
            human_cell_donor_or_batch=donor,target=target,
            condition=condition,technical_set=set_id,
            normalized_targeted_luminal_signal=targeted,
            normalized_size_matched_nontarget_signal=nonspecific,
            normalized_competition_blocked_luminal_signal=blocked,
            targeted_minus_nonbinding_signal=net,
            fold_over_nonbinding=(targeted/nonspecific if nonspecific>0 else ""),
            competition_reduction_fraction=(competitor_reduction if
                                            competitor_reduction is not None else ""),
            maximal_dextran_permeability=max(x["permeability"] for x in arms.values()),
            barrier_acceptable="PASS_LOCKED_THRESHOLD",
            human_ischemic_brain_vessel_binding="NOT_TESTED",
            PET_specificity="NOT_TESTED"))
    by_donor=defaultdict(list)
    for r in experimental:
        by_donor[(r["human_cell_donor_or_batch"],r["target"],r["condition"])].append(r)
    donor_results=[]
    for (donor,target,condition),items in sorted(by_donor.items()):
        donor_results.append(dict(
            donor_batch_id=donor,target=target,condition=condition,
            technical_triples=len(items),
            median_normalized_target_minus_nonbinding_signal=
                statistics.median(x["targeted_minus_nonbinding_signal"] for x in items),
            median_competitor_reduction_fraction=statistics.median(
                x["competition_reduction_fraction"] for x in items
                if x["competition_reduction_fraction"]!="") if any(
                x["competition_reduction_fraction"]!="" for x in items) else "",
            independent_biological_replicate=1,
            assay_scope="HUMAN_BBB_IN_VITRO_NOT_HUMAN_ISCHEMIC_STROKE"))
    summary=[]
    for target in sorted(TARGETS):
        for condition in sorted(CONDITIONS):
            g=[r for r in donor_results if r["target"]==target and
               r["condition"]==condition]
            summary.append(dict(target=target,condition=condition,
                n_independent_human_cell_donors_or_batches=len(g),
                normalized_binding_difference_median_across_donors=(
                    statistics.median(x["median_normalized_target_minus_nonbinding_signal"]
                                      for x in g) if g else ""),
                competition_reduction_median_across_donors=(
                    statistics.median(x["median_competitor_reduction_fraction"] for x in g
                                      if x["median_competitor_reduction_fraction"]!="")
                    if any(x["median_competitor_reduction_fraction"]!="" for x in g) else ""),
                status=("INSUFFICIENT_BIOLOGICAL_REPLICATION_DESCRIPTIVE" if len(g)<3
                        else "IN_VITRO_FEASIBILITY_ONLY_NOT_STROKE_PROTEIN_OR_PET")))
    destination.mkdir(parents=True,exist_ok=True)
    write(destination/"24_luminal_binding_matched_chip_triples.csv",experimental)
    write(destination/"24_human_BBB_donor_batch_level_binding_summaries.csv",donor_results)
    write(destination/"24_live_luminal_binding_feasibility_overview.csv",summary)
    metadata={
      "n_assay_chip_records":len(rows),"n_matched_ligand_triples":len(experimental),
      "n_biological_cell_donor_target_condition_groups":len(donor_results),
      "n_technical_paired_chips_per_set":3,
      "model":"Live nonpermeabilized perfused BBB chip, NOT tissue taken from human ischemic stroke",
      "controls":"size/exposure matched nonbinding and competition arms, locked barrier permeability and QC",
      "lumen_binding_not_equivalent_to_human_ischemic_brain_accessibility":True,
      "clinical_protein_and_pet_validated":False,
      "bio_unit":"independent cell donor/differentiation batch",
      "chip_unit":"technical",
      "synthetic_data_may_never_be_interpreted_as_experiment":True
    }
    (destination/"24_live_BBB_luminal_binding_analysis_manifest.json").write_text(
        json.dumps(metadata,indent=2)+"\n",encoding="utf-8")
    return metadata

def selftest():
    # Synthetic inputs deleted after smoke-test, NEVER archived as biological evidence.
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);infile=root/"chip.csv";controls=root/"locked.csv"
        write(controls,[dict(label_batch="L1",
            max_acceptable_dextran_70kDa_permeability=.001,
            barrier_limit_locked_before_condition_unblinding="True",
            positive_live_surface_control_pass="True")])
        fake=[]
        for donor in range(1,4):
            for target in sorted(TARGETS):
                for arm,signal in (("target_ligand",100),
                                   ("size_matched_nonbinder",20),
                                   ("competitor_blocked_target",30)):
                    fake.append(dict(
                       donor_batch_id=f"CB{donor}",technical_set_id="set1",
                       chip_id=f"C{donor}_{target}_{arm}",target=target,
                       condition="TNF_activation",assay_arm=arm,assay="LIVE_LUMINAL_FLOW",
                       label_batch="L1",ligand_size_kDa=55,input_probe_AUC=100,
                       luminal_retained_probe_signal=signal,
                       luminal_endothelial_area=10,extravasated_probe_signal=1,
                       permeability_70kDa=.0002,
                       viability_pass="True",barrier_structure_pass="True",
                       probe_added_from_lumen="True",nonpermeabilized="True",
                       images_blinded="True"))
        write(infile,fake)
        report=run(infile,controls,root/"okay")
        require(report["n_matched_ligand_triples"]==6,
                "Treated technical assays as independent biological chips")
        invalid=[dict(x) for x in fake]
        invalid[-1]["permeability_70kDa"]=.1
        write(infile,invalid)
        try: run(infile,controls,root/"no_leak")
        except ValueError:pass
        else:raise AssertionError("Leaky BBB passed live binding gate")
        invalid=[dict(x) for x in fake]
        invalid.pop()
        write(infile,invalid)
        try: run(infile,controls,root/"no_blocker")
        except ValueError:pass
        else:raise AssertionError("Missing blocker passed gate")
        print("LIVE_HUMAN_BBB_BINDER_ANALYSIS_SELFTEST PASS",
              "synthetic only; donor-level distinct n=3; reject leak/missing blocker")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--chips",type=Path)
    p.add_argument("--locked-barrier-controls",type=Path)
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    if args.self_test:
        selftest()
    else:
        require(args.chips and args.locked_barrier_controls and args.out,
                "Supply --chips, --locked-barrier-controls, --out for real data")
        report=run(args.chips,args.locked_barrier_controls,args.out)
        print("LIVE_HUMAN_BBB_BINDER_AUDIT PASS",json.dumps(report),flush=True)
