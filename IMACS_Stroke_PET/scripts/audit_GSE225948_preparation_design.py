#!/usr/bin/env python3
"""Independent source-aware GSE225948 design, composition and stability audit.

Inputs: immutable archived preparation-level results and the PASS artifact from
the separately executed full-matrix audit. This script does not use cells as
replicates and never labels pooled preparations as independent animals.

Caution: label-permutation p-values are arithmetically reproduced, not endorsed
as calibrated disease p-values when processing batches are not exchangeable.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "preparation_design_qc"

GROUPS = {
    "Sham": [f"GSM70608{i:02d}" for i in range(15, 19)],
    "Day2": [f"GSM70608{i:02d}" for i in range(19, 23)],
    "Day14": [f"GSM70608{i:02d}" for i in range(23, 27)],
}
ALL_IDS = [x for xs in GROUPS.values() for x in xs]


def rows(path: Path) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def put_csv(path: Path, data: list[dict]):
    if not data:
        raise AssertionError(f"No rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def close(a, b, tol=1e-9):
    return math.isclose(float(a), float(b), rel_tol=tol, abs_tol=tol)


def fail(condition, message):
    if not condition:
        raise AssertionError(message)


def group_of(gsm):
    return next((g for g, ids in GROUPS.items() if gsm in ids), None)


def keyed(data, keys):
    result = {}
    for row in data:
        key = row[keys[0]] if len(keys) == 1 else tuple(row[k] for k in keys)
        fail(key not in result, f"Duplicate row for {key}")
        result[key] = row
    return result


def permutation(vals):
    fail(len(vals) == 8 and all(math.isfinite(x) for x in vals), "Invalid balanced comparison")
    observed = sum(vals[4:]) / 4 - sum(vals[:4]) / 4
    ps = []
    for chosen in itertools.combinations(range(8), 4):
        chosen = set(chosen)
        case = sum(vals[i] for i in chosen) / 4
        ctrl = sum(vals[i] for i in range(8) if i not in chosen) / 4
        ps.append(case - ctrl)
    fail(len(ps) == 70, "Expected 70 assignments")
    p = sum(abs(x) >= abs(observed) - 1e-12 for x in ps) / len(ps)
    return observed, p


def bh(p_values):
    ordered = sorted(range(len(p_values)), key=lambda i: p_values[i])
    out = [None] * len(ordered)
    lowest = 1.0
    for rank in range(len(ordered), 0, -1):
        idx = ordered[rank - 1]
        lowest = min(lowest, p_values[idx] * len(ordered) / rank)
        out[idx] = lowest
    return out


def main(qc_dir):
    OUT.mkdir(parents=True, exist_ok=True)
    day = RES / "independent_day2"
    late = RES / "independent_late"
    manifests = [rows(day / "11_inputs_checksums.csv"),
                 rows(late / "12_inputs_checksums.csv")]
    source = {}
    for manifest in manifests:
        for r in manifest:
            gsm, typ = r["gsm"], r["type"]
            fail(gsm in ALL_IDS and typ in {"counts", "metadata"}, "Unrecognized source")
            fail(re.fullmatch(r"[a-f0-9]{64}", r["sha256"]) is not None, "Invalid SHA256")
            fail(r["file"].startswith(gsm + "_Brain_GR"), "Unexpected file/prefix")
            fail("ftp.ncbi.nlm.nih.gov/geo/samples/" in r["url"], "Unexpected GEO host")
            key = gsm, typ
            if key in source:
                fail((source[key]["sha256"], source[key]["file"], source[key]["url"]) ==
                     (r["sha256"], r["file"], r["url"]), f"Source discrepancy {key}")
            else:
                source[key] = r
    fail(len(source) == 24, "Expected 12 distinct matrix/metadata pairs")
    for group, ids in GROUPS.items():
        for gsm in ids:
            fail((gsm, "counts") in source and (gsm, "metadata") in source,
                 f"Missing source {gsm}")

    full_qc = keyed(rows(qc_dir / "13_GSE225948_full_matrix_normalization_audit.csv"), ["gsm"])
    fail(set(full_qc) == set(ALL_IDS), "Missing or extra normalization result")
    for gsm, qc in full_qc.items():
        fail(qc["status"] == "PASS", f"Normalization gate fails for {gsm}")
        fail(close(qc["ec_median_sum_over_nCount"], 1.0, 1e-8), f"nCount mismatch {gsm}")
        fail(close(qc["ec_fraction_within_1pct"], 1.0), f"EC matrix totals mismatch {gsm}")
        fail(close(qc["ec_nFeature_exact_fraction"], 1.0), f"nFeature mismatch {gsm}")
        fail(int(qc["genes"]) == 19704, "Unanticipated feature dimension")

    comp_a = rows(day / "11_author_annotated_population_counts.csv")
    comp_b = rows(late / "12_author_annotated_population_counts.csv")
    comp = defaultdict(dict)
    for records in (comp_a, comp_b):
        current = defaultdict(dict)
        for r in records:
            gsm, label, n = r["gsm"], r["author_parent"], int(r["n_cells"])
            fail(group_of(gsm) is not None and n >= 0, "Invalid cell annotation")
            fail(label not in current[gsm], "Duplicate lineage row")
            current[gsm][label] = n
        for gsm, lineages in current.items():
            if gsm in comp:
                fail(comp[gsm] == lineages, f"Cross-notebook lineage mismatch {gsm}")
            comp[gsm] = lineages
    fail(set(comp) == set(ALL_IDS), "Missing lineage-annotation group")

    cell_qc = {}
    for p, prefix in ((day, "11"), (late, "12")):
        for r in rows(p / f"{prefix}_matrix_format_and_barcode_QC.csv"):
            gsm = r["gsm"]
            fail(int(r["unmatched"]) == 0 and
                 int(r["cells_expression"]) == int(r["matched_metadata"]),
                 f"Cell barcode mismatch {gsm}")
            fail(int(r["duplicate_feature_names"]) == 0, f"Duplicate marker features {gsm}")
            fail(int(r["cells_expression"]) == sum(comp[gsm].values()),
                 f"Lineage totals do not sum to library {gsm}")
            if gsm in cell_qc:
                fail(int(cell_qc[gsm]["cells_expression"]) ==
                     int(r["cells_expression"]), "Cross-run cell total disagreement")
            cell_qc[gsm] = r
    fail(set(cell_qc) == set(ALL_IDS), "Cell QC incomplete")

    whole = {}
    by_subtype = {}
    actual_genes = None
    for p, prefix, ids, group_a, group_b in (
        (day, "11", GROUPS["Sham"] + GROUPS["Day2"], "Sham", "Stroke"),
        (late, "12", GROUPS["Day2"] + GROUPS["Day14"], "Day2", "Day14"),
    ):
        records = rows(p / f"{prefix}_day2_EC_gene_by_preparation.csv")
        matched = [r for r in records if r["status"] == "TESTABLE"]
        genes = sorted({r["gene"] for r in matched})
        fail(len(genes) == 28, "Expected 28 testable frozen-panel genes")
        fail(actual_genes is None or genes == actual_genes, "Panel changed across contrasts")
        actual_genes = genes
        v = keyed(matched, ["gsm", "gene"])
        fail(len(v) == 8 * len(genes), "Missing a preparation/gene")
        for gsm in ids:
            fail(int(next(v[gsm, gene] for gene in genes)["n_author_EC_cells"]) ==
                 comp[gsm]["EC"], f"EC count mismatch {gsm}")
            for gene in genes:
                item = v[gsm, gene]
                fail(item["condition"] == (group_a if gsm in ids[:4] else group_b),
                     f"Incorrect treatment labels {gsm}")
                fail(int(item["n_author_EC_cells"]) == comp[gsm]["EC"], "EC cell inconsistency")
                val = float(item["mean_cell_log1p_10000"])
                fail(math.isfinite(val) and val >= -1e-10, "Invalid panel expression")
                if gsm in GROUPS["Day2"]:
                    k = gsm, gene
                    if k in whole:
                        fail(close(whole[k], val, 1e-10), f"Day2 independent-run mismatch {k}")
                    else:
                        whole[k] = val

        observed = keyed(rows(p / f"{prefix}_GSE225948_day2_EC_exact_permutation.csv"),
                         ["gene"])
        fail(set(observed) == set(genes), "Saved permutation panel mismatch")
        ec1_records = [r for r in rows(p / f"{prefix}_author_EC1_EC2_sample_gene_descriptive.csv")
                       if r["gene"] in genes]
        subtype_by_key = keyed(ec1_records, ["gsm", "author_subtype", "gene"])
        ec1 = keyed([r for r in rows(p / f"{prefix}_within_author_EC_subtype_sensitivity.csv")
                     if r["author_subtype"] == "EC1"], ["gene"])
        fail(set(ec1) == set(genes), "Missing EC1 sensitivity results")

        subtype_counts = {}
        for gsm in ids:
            n = {}
            for subtype in ("EC1", "EC2"):
                counts = {int(subtype_by_key[gsm, subtype, g]["n_cells"])
                          for g in genes if (gsm, subtype, g) in subtype_by_key}
                fail(len(counts) == 1, f"Inconsistent {subtype} sample counts {gsm}")
                n[subtype] = next(iter(counts))
            fail(sum(n.values()) <= comp[gsm]["EC"], f"Subtype total exceeds EC {gsm}")
            subtype_counts[gsm] = n

        group_rows = []
        p_list = []
        recomputed = []
        for gene in genes:
            vals = [float(v[gsm, gene]["mean_cell_log1p_10000"]) for gsm in ids]
            delta, p_value = permutation(vals)
            truth = observed[gene]
            delta_col = ("stroke_minus_sham_mean_log1p_10000" if prefix == "11"
                         else "day14_minus_day2_mean_log1p_10000")
            fail(close(delta, truth[delta_col], 1e-10) and
                 close(p_value, truth["exact_two_sided_p"], 1e-10),
                 f"Recomputed original test differs for {gene} ({prefix})")
            p_list.append(p_value)
            recomputed.append((gene, vals, delta, p_value))
        recalculated_q = bh(p_list)
        for (gene, vals, delta, p_value), q in zip(recomputed, recalculated_q):
            fail(close(q, observed[gene]["BH_q_across_frozen_panel"], 1e-10),
                 f"Frozen panel BH mismatch {prefix} {gene}")
            leave = []
            for i in range(4):
                for j in range(4, 8):
                    a = [vals[k] for k in range(4) if k != i]
                    b = [vals[k] for k in range(4, 8) if k != j]
                    leave.append(sum(b) / 3 - sum(a) / 3)
            ec1_vals = [float(subtype_by_key[gsm, "EC1", gene]["mean_log1p_10000"])
                        for gsm in ids]
            ec1_delta, ec1_p = permutation(ec1_vals)
            ec1_saved = ec1[gene]
            fail(ec1_saved["status"] == "TESTABLE_4v4_pooled_preps", "EC1 test not available")
            fail(close(ec1_delta, ec1_saved["delta_aggregate_mean_log1p_10000"], 1e-10)
                 and close(ec1_p, ec1_saved["exact_two_sided_p"], 1e-10),
                 f"EC1 reference mismatch for {prefix}/{gene}")
            flips = sum(x * delta < -1e-12 for x in leave)
            group_rows.append({
                "contrast": "StrokeDay2_vs_ShamDay2" if prefix == "11" else "StrokeDay14_vs_StrokeDay2",
                "gene": gene,
                "whole_ec_delta_mean_log1p10000": delta,
                "archived_exact_p_conditional_on_exchangeability": p_value,
                "archived_BH_q_frozen_panel": q,
                "ec1_delta_mean_log1p10000": ec1_delta,
                "ec1_archived_exact_p_conditional_on_exchangeability": ec1_p,
                "ec1_BH_q_frozen_panel": ec1_saved["BH_q_within_subtype_frozen_panel"],
                "same_direction_whole_ec_and_ec1": (delta * ec1_delta > 0),
                "leave_one_per_arm_min_delta": min(leave),
                "leave_one_per_arm_max_delta": max(leave),
                "leave_one_per_arm_direction_reversals_of_16": flips,
                "leave_one_per_arm_n_variants": 16,
                "interpretation": "DESCRIPTIVE_ONLY_BATCH_EXCHANGEABILITY_NOT_VALIDATED",
            })
        by_subtype[prefix] = subtype_counts
        # Sample-specific capture is NOT tissue abundance.
        by_subtype[prefix + "_results"] = group_rows

    output = []
    for gsm in ALL_IDS:
        group = group_of(gsm)
        n_ec = comp[gsm]["EC"]
        parts = by_subtype["11" if gsm in GROUPS["Sham"] else "12"][gsm]
        label = re.search(r"_Brain_(GR\d{6})_", source[gsm, "counts"]["file"])
        fail(label is not None, f"Missing preparation key {gsm}")
        output.append({
            "gsm": gsm, "experimental_group": group,
            "original_geo_series": "GSE225948",
            "n_deposited_preparations_per_group": 4,
            "n_individual_animals_traced": "UNRESOLVED",
            "source_reports_hemispheres_per_pool": "4-5",
            "preparation_identifier_date_like_not_verified_batch": label.group(1),
            "all_captured_author_labeled_cells": sum(comp[gsm].values()),
            "author_ec_cells": n_ec,
            "ec_fraction_of_sorted_captured_cells": n_ec / sum(comp[gsm].values()),
            "author_ec1_cells": parts["EC1"],
            "author_ec2_cells": parts["EC2"],
            "ec1_fraction_within_captured_ec": parts["EC1"] / n_ec,
            "ec2_fraction_within_captured_ec": parts["EC2"] / n_ec,
            "other_ec_subtype_fraction": (n_ec - sum(parts.values())) / n_ec,
            "full_matrix_normalization": full_qc[gsm]["status"],
            "nCount_median_matrix_ratio": full_qc[gsm]["ec_median_sum_over_nCount"],
        })

    stability = by_subtype["11_results"] + by_subtype["12_results"]
    fail(len(stability) == 56, "Expected 2 contrasts x 28 testable genes")
    report = {}
    for group, ids in GROUPS.items():
        r = [x for x in output if x["gsm"] in ids]
        report[group] = {
            "preparations": len(r),
            "captured_ec_cells_total_not_biological_n": sum(x["author_ec_cells"] for x in r),
            "captured_ec_fraction_median": sorted(x["ec_fraction_of_sorted_captured_cells"] for x in r)[len(r)//2 - 1:len(r)//2 + 1],
            "ec2_fraction_range": [min(x["ec2_fraction_within_captured_ec"] for x in r),
                                   max(x["ec2_fraction_within_captured_ec"] for x in r)],
            "ec2_fraction_mean_descriptive": sum(x["ec2_fraction_within_captured_ec"] for x in r) / len(r),
            "date_like_preparation_labels": [x["preparation_identifier_date_like_not_verified_batch"] for x in r],
        }
    summaries = {}
    for code, name in (("11", "StrokeDay2_vs_ShamDay2"), ("12", "StrokeDay14_vs_StrokeDay2")):
        r = by_subtype[code + "_results"]
        summaries[name] = {
            "testable_frozen_genes": len(r),
            "BH_q_below_0_05": sum(float(x["archived_BH_q_frozen_panel"]) < .05 for x in r),
            "ec1_BH_q_below_0_05": sum(float(x["ec1_BH_q_frozen_panel"]) < .05 for x in r),
            "contrasts_with_leave_one_per_arm_direction_reversal":
                sum(x["leave_one_per_arm_direction_reversals_of_16"] > 0 for x in r),
            "sign_disagreements_EC1_vs_whole":
                [x["gene"] for x in r if x["whole_ec_delta_mean_log1p10000"] *
                 x["ec1_delta_mean_log1p10000"] < -1e-12],
        }
    result = {
        "dataset": "GSE225948", "qc_source": "executed GitHub Actions run 38427680867",
        "source_files_sha256_compatible": True, "barcode_and_annotation_totals_match": True,
        "full_matrix_nCount_and_nFeature_agreement": "PASS (12/12 preparations)",
        "source_accessions": ALL_IDS,
        "design": {
            "unit": "pooled sorted preparation, 4-5 hemispheres per deposited experiment (GEO)",
            "single_animal_to_cell_mapping": "NOT AVAILABLE IN INPUT MANIFEST; individual-animal DE impossible",
            "sham_day14": "ABSENT; day14 vs day2 is stroke-only comparison",
            "batch_blocking": "No pair of groups shares an identical preparation ID; processing batch overlap is not established",
            "date_like_identifiers": "Filename GRYYMMDD-like strings; temporal interpretation NOT independently verified",
            "group_exchangeability": "UNVERIFIED; archival exact permutation p and BH are conditional arithmetic only",
            "annotation": "Author-provided parent==EC, EC1 and EC2. Independent lineage validation NOT performed",
            "capture": "FACS-enriched cell samples, cell fractions are NOT tissue cell abundance",
            "normalization": "Fractional non-integer released counts. No DESeq2 on these counts",
            "pet": "RNA results do not establish protein abundance, vascular luminal access, binding, or PET feasibility",
        },
        "preparation_group_summary": report,
        "effect_summary": summaries,
        "method": "Strict provenance and reference p/BH recomputation; 16 deterministic leave-one-preparation-per-arm descriptive contrasts for each gene",
    }
    put_csv(OUT / "14_GSE225948_preparation_and_subtype_design.csv", output)
    put_csv(OUT / "14_GSE225948_gene_leave_one_preparation_sensitivity.csv", stability)
    (OUT / "14_GSE225948_design_and_stability_audit.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print("DESIGN_AUDIT PASS", json.dumps(summaries), flush=True)
    print("GSE225948 sample capture:", json.dumps(report), flush=True)
    print("Scientific gate: RESTRICTED, NOT biological target validation", flush=True)


if __name__ == "__main__":
    argp = argparse.ArgumentParser()
    argp.add_argument("--qc-dir", type=Path, required=True)
    args = argp.parse_args()
    main(args.qc_dir)
