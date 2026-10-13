#!/usr/bin/env python3
"""Audit rat GSE250245 pooled preparation results without pseudoreplication.

Reassesses archived comparison eligibility. A sampled population with >=50
cells cannot rescue a gene missing from deposited features. No group-level
differential-expression p-values are calculated from pooled matrices.
"""
from __future__ import annotations
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "results" / "rat"
OUT = BASE / "results" / "rat_provenance_qc"
OUT.mkdir(parents=True, exist_ok=True)


def read(path):
    with (DATA / path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def require(x, msg):
    if not x:
        raise AssertionError(msg)


def write(name, records):
    require(records, "Empty output " + name)
    with (OUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    source = read("03_rat_GEO_design_manifest.csv")
    required = [
        ("GSM7976206", "Sham_L", "Sham", "contralateral_or_left", 4),
        ("GSM7976207", "Sham_R", "Sham", "ipsilateral_or_right", 4),
        ("GSM7976208", "mMCAO_contra", "moderate", "contralateral", 3),
        ("GSM7976209", "mMCAO_ipsi", "moderate", "ipsilateral", 3),
        ("GSM7976210", "sMCAO_contra", "severe", "contralateral", 4),
        ("GSM7976211", "sMCAO_ipsi", "severe", "ipsilateral", 4),
    ]
    require(len(source) == 6, "Expected exactly 6 GEO group-hemisphere libraries")
    by_id = {r["gsm"]: r for r in source}
    require(len(by_id) == 6, "Duplicate or unexpected accession")
    for gsm, label, severity, hemisphere, n in required:
        r = by_id[gsm]
        require((r["prefix"], r["severity"], r["hemisphere"],
                 int(r["reported_group_n_rats"])) ==
                (label, severity, hemisphere, n),
                "GEO design inconsistency for " + gsm)

    checksum = read("03_rat_SHA256.csv")
    require(len(checksum) == 18, "Expected three original files per GEO sample")
    files = set()
    for r in checksum:
        require(r["gsm"] in by_id and r["file"].startswith(r["gsm"] + "_"),
                "Invalid source file")
        require(int(r["bytes"]) > 0 and len(r["sha256"]) == 64 and
                all(ch in "0123456789abcdef" for ch in r["sha256"]),
                "Corrupt source provenance")
        files.add((r["gsm"], r["file"]))
    require(len(files) == 18, "Duplicated source file provenance")

    qc = {r["gsm"]: r for r in read("03_rat_48h_QC_audit.csv")}
    require(set(qc) == set(by_id), "QC does not cover all six pooled matrices")
    sample_rows = []
    for gsm, name, sev, hemi, n in required:
        r = qc[gsm]
        before, after = int(r["n_cells_before"]), int(r["n_cells_after"])
        require(before >= after >= 1, "QC negative or impossible")
        sample_rows.append({
            "gsm": gsm, "source_group": name, "rat_group_size_reported": n,
            "individual_rat_to_barcode_map": "UNAVAILABLE_IN_ARCHIVED_INPUTS",
            "nuclei_before_inhouse_qc": before, "nuclei_after_inhouse_qc": after,
            "nuclei_kept_fraction": after / before,
            "median_features_after_inhouse_qc": r["median_genes_after"],
            "median_mt_fraction_after_inhouse_qc": r["median_pct_mt_after"],
        })

    missing = read("03_rat_missing_gene_symbol_audit.csv")
    absent = defaultdict(set)
    for r in missing:
        require(r["gsm"] in by_id and r["not_found"] == "True",
                "Unexpected feature diagnostic")
        absent[r["gsm"]].add(r["gene"])
    require(set(absent) == set(by_id), "Missing gene audit incomplete")
    require(len(set(tuple(sorted(v)) for v in absent.values())) == 1,
            "Rat gene feature sets disagree on frozen-panel missing genes")
    absent_genes = sorted(next(iter(absent.values())))

    expression = read("03_rat_48h_candidate_expression_by_group_and_hemisphere.csv")
    expr = {}
    populations = {}
    for r in expression:
        key = (r["prefix"], r["celltype"], r["gene"])
        require(key not in expr, "Duplicate celltype/gene expression row")
        require(r["gsm"] in by_id and by_id[r["gsm"]]["prefix"] == r["prefix"],
                "Incorrect sample-to-condition mapping")
        expr[key] = r
        pop = (r["gsm"], r["celltype"])
        n = int(r["n_cells"])
        require(n >= 1, "Empty recorded population")
        if pop in populations:
            require(populations[pop] == n, "Population count differs across genes")
        populations[pop] = n
        if r["gene"] in absent[r["gsm"]]:
            require(r["mean_log1p"] == "" and r["detection_fraction"] == "",
                    "Absent feature falsely represented as observed expression")
        else:
            y, detect = float(r["mean_log1p"]), float(r["detection_fraction"])
            require(math.isfinite(y) and y >= -1e-9 and
                    math.isfinite(detect) and 0 <= detect <= 1,
                    "Invalid measured gene expression")
    for gsm, name, *_ in required:
        require(sum(c for (key, typ), c in populations.items() if key == gsm) ==
                int(qc[gsm]["n_cells_after"]),
                "Annotated nuclei fail to sum to post-QC sample count")

    revised = []
    old_false_pass = 0
    low_endothelial = []
    original_comp = read("03_rat_group_comparisons_descriptive_NO_PVALUES.csv")
    require(all("p_value" not in r and "q_value" not in r for r in original_comp),
            "Unexpected inferential column for pooled rat data")
    for r in original_comp:
        gene = r["gene"]
        label = r["celltype"]
        case, control = r["case"], r["control"]
        require((case, label, gene) in expr and (control, label, gene) in expr,
                "Comparison not found in archived source")
        ca, co = expr[(case, label, gene)], expr[(control, label, gene)]
        n_case, n_ctrl = int(r["case_n_cells"]), int(r["control_n_cells"])
        require(n_case == int(ca["n_cells"]) and n_ctrl == int(co["n_cells"]),
                "Comparison cell count mismatch")
        absent_on_source = ca["mean_log1p"] == "" or co["mean_log1p"] == ""
        if absent_on_source:
            gate = "BLOCK_MISSING_FEATURE"
            require(r["log_mean_difference"] == "" and
                    r["detection_fraction_difference"] == "",
                    "Missing feature has numerical comparison")
            if r["comparison_qc"] == "PASS_COUNTS_ONLY_DESCRIPTIVE":
                old_false_pass += 1
        elif min(n_case, n_ctrl) < 50:
            gate = "BLOCK_LOW_CAPTURE"
        else:
            gate = "DESCRIPTIVE_ONLY_NO_RAT_LEVEL_REPLICATION"
        if not absent_on_source:
            delta = float(ca["mean_log1p"]) - float(co["mean_log1p"])
            require(math.isclose(delta, float(r["log_mean_difference"]),
                                 abs_tol=2e-6, rel_tol=2e-6),
                    "Archived expression difference does not reproduce")
        revised.append({
            "celltype": label, "gene": gene, "case": case, "control": control,
            "case_cells_post_QC": n_case, "control_cells_post_QC": n_ctrl,
            "source_gene_feature_present": not absent_on_source,
            "revised_scientific_gate": gate,
            "descriptive_log_mean_difference": r["log_mean_difference"]
                if not absent_on_source else "",
            "original_comparison_qc": r["comparison_qc"],
            "inference": "NO_RAT_LEVEL_PVALUE_OR_PET_CLAIM",
        })

    rat_cells = sum(int(r["n_cells_after"]) for r in qc.values())
    for gsm, name, *_ in required:
        k = (gsm, "Endothelial")
        n = populations.get(k, 0)
        if n < 50:
            low_endothelial.append({"gsm": gsm, "source_group": name, "EC_nuclei": n})
        sample_rows[[x["gsm"] for x in sample_rows].index(gsm)]["inhouse_marker_endothelial_nuclei"] = n

    result = {
        "dataset": "GSE250245", "organism": "Rattus norvegicus",
        "deposited_matrices": 6,
        "source_groups": "4 sham rats, 3 moderate MCAO rats, 4 severe MCAO rats, bilateral hemispheres",
        "specimen_to_rat_mapping_available_for_inhouse_pseudobulk": False,
        "independent_rat_level_differential_expression_possible_from_archived_matrices": False,
        "rat_study_design_authority": "GEO GSE250245 / PMID 39043661",
        "inhouse_qc_retained_nuclei_all_matrices": rat_cells,
        "author_publication_clustered_nuclei": 68616,
        "difference_inhouse_minus_author": rat_cells - 68616,
        "difference_interpretation": "Distinct preprocessing/cell-filtering and annotation; requires recovery of original QC/metadata before asserting lineage concordance",
        "inhouse_annotation": "custom broad marker-score classifier, not independently validated",
        "source_missing_frozen_genes_all_six_groups": absent_genes,
        "missing_gene_vs_biological_nonexpression": "NOT EQUIVALENT; no fuzzy substitution",
        "source_comparison_rows": len(original_comp),
        "original_count_only_PASS_rows_with_missing_gene": old_false_pass,
        "comparison_gate_counts": dict(sorted(__import__("collections").Counter(
            r["revised_scientific_gate"] for r in revised).items())),
        "low_endothelial_capture_lt50": low_endothelial,
        "interpretation": "Exploratory pooled cross-model context only; cannot validate vascular stroke targets or PET suitability",
    }
    write("15_rat_source_population_QC.csv", sample_rows)
    write("15_rat_descriptive_comparisons_reclassified.csv", revised)
    (OUT / "15_rat_source_annotation_and_missing_gene_audit.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("RAT_AUDIT PASS", json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
