#!/usr/bin/env python3
"""Audit full GSE225948 fractional matrices against author nCount_RNA metadata.

This is a source-integrity/normalization gate, NOT a differential-expression test.
Never interpret a discrepancy as biological absence or use cell counts as replicates.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / "source_normalization_qc"
RAW = Path(os.environ.get("IMACS_QC_RAW", "/tmp/imacs_gse225948_qc"))
RAW.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

MANIFESTS = [
    BASE / "results" / "independent_day2" / "11_inputs_checksums.csv",
    BASE / "results" / "independent_late" / "12_inputs_checksums.csv",
]
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=Retry(
    total=4, connect=4, read=4, backoff_factor=2,
    status_forcelist=[429, 500, 502, 503, 504]
)))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_manifest() -> pd.DataFrame:
    rows = pd.concat([pd.read_csv(p) for p in MANIFESTS], ignore_index=True)
    rows = rows[rows["type"].isin(["counts", "metadata"])].copy()
    assert rows["gsm"].nunique() == 12, "Expected exactly 12 GEO preparations"
    for (gsm, typ), group in rows.groupby(["gsm", "type"]):
        assert len(group["sha256"].unique()) == 1, (gsm, typ, "checksum conflict")
        assert len(group["url"].unique()) == 1, (gsm, typ, "URL conflict")
    rows = rows.drop_duplicates(["gsm", "type"]).sort_values(["gsm", "type"])
    assert len(rows) == 24
    for row in rows.itertuples(index=False):
        parsed = urlparse(row.url)
        assert parsed.scheme == "https" and parsed.hostname == "ftp.ncbi.nlm.nih.gov"
        assert f"/{row.gsm}/suppl/{row.gsm}_" in parsed.path
        assert str(row.sha256) and len(str(row.sha256)) == 64
    return rows


def get_file(row) -> Path:
    target = RAW / str(row.file)
    expected = str(row.sha256)
    if target.exists() and sha256(target) == expected:
        return target
    if target.exists():
        target.unlink()
    temp = target.with_suffix(target.suffix + ".partial")
    try:
        with session.get(row.url, timeout=(30, 300), stream=True) as response:
            response.raise_for_status()
            with temp.open("wb") as stream:
                for block in response.iter_content(1 << 20):
                    if block:
                        stream.write(block)
        assert sha256(temp) == expected, f"SHA256 mismatch: {row.gsm} {row.type}"
        temp.replace(target)
    finally:
        if temp.exists():
            temp.unlink()
    return target


def compare_sums(gsm: str, count_path: Path, meta_path: Path) -> dict:
    meta = pd.read_csv(meta_path, compression="gzip", index_col=0)
    assert meta.index.is_unique, (gsm, "duplicate metadata barcodes")
    assert {"nCount_RNA", "nFeature_RNA", "parent", "sub.celltype"}.issubset(meta.columns)
    meta.index = meta.index.astype(str)
    lib = pd.to_numeric(meta["nCount_RNA"], errors="raise")
    features = pd.to_numeric(meta["nFeature_RNA"], errors="raise")
    assert np.isfinite(lib).all() and (lib > 0).all()
    assert np.isfinite(features).all() and (features >= 0).all()

    sums = np.zeros(len(meta), dtype=np.float64)
    detected = np.zeros(len(meta), dtype=np.int64)
    seen = set()
    total_genes = 0
    noninteger_seen = False
    ordered_cols = None
    for chunk in pd.read_csv(count_path, compression="gzip", index_col=0, chunksize=128):
        cols = chunk.columns.astype(str)
        if ordered_cols is None:
            ordered_cols = cols
            assert len(cols) == len(meta) and cols.is_unique
            assert set(cols) == set(meta.index), (gsm, "count/metadata barcode mismatch")
            order = pd.Index(cols).get_indexer(meta.index)
            assert (order >= 0).all()
        else:
            assert ordered_cols.equals(cols), (gsm, "column drift")
        gene_names = chunk.index.astype(str)
        assert not gene_names.has_duplicates
        assert not seen.intersection(gene_names)
        seen.update(gene_names)
        vals = chunk.to_numpy(dtype=np.float64, copy=False)
        assert np.isfinite(vals).all() and (vals >= 0).all(), (gsm, "invalid matrix values")
        noninteger_seen |= bool(np.any(np.abs(vals - np.rint(vals)) > 1e-7))
        vals = vals[:, order]
        sums += vals.sum(axis=0, dtype=np.float64)
        detected += np.count_nonzero(vals > 0, axis=0)
        total_genes += len(chunk)

    assert total_genes > 1000 and ordered_cols is not None
    ratio = sums / lib.to_numpy(float)
    fdelta = detected - features.to_numpy(float)
    ec = meta["parent"].astype(str).eq("EC").to_numpy()
    assert ec.sum() >= 10, (gsm, "too few annotated endothelial cells")
    if not (np.isfinite(ratio).all() and (sums > 0).all()):
        raise AssertionError((gsm, "nonfinite or zero source library sums"))
    def summary(mask):
        x = ratio[mask]
        return {
            "n_cells": int(mask.sum()),
            "median_sum_over_nCount": float(np.median(x)),
            "p05_sum_over_nCount": float(np.percentile(x, 5)),
            "p95_sum_over_nCount": float(np.percentile(x, 95)),
            "fraction_within_1pct": float(np.mean(np.abs(x - 1) <= .01)),
            "fraction_within_5pct": float(np.mean(np.abs(x - 1) <= .05)),
        }
    ec_summary = summary(ec)
    all_summary = summary(np.ones(len(meta), dtype=bool))
    # Conservative gate: denominators need to represent this matrix's library totals
    # to support the current per-cell log1p(1e4*x/nCount_RNA) interpretation.
    status = "PASS" if ec_summary["fraction_within_5pct"] >= .95 else "REVIEW"
    return {
        "gsm": gsm, "status": status, "genes": total_genes,
        "fractional_values_present": bool(noninteger_seen),
        "ec": ec_summary, "all_cells": all_summary,
        "ec_nFeature_exact_fraction": float(np.mean(fdelta[ec] == 0)),
        "ec_nFeature_median_difference": float(np.median(fdelta[ec])),
        "interpretation": (
            "PASS means numerical library denominator agreement only, not biological validation."
            if status == "PASS" else
            "Matrix column sums disagree with nCount_RNA for >5% of author EC cells; "
            "revisit normalization and all previously derived expression effects."
        ),
    }


def main():
    manifest = load_manifest()
    results = []
    for gsm, group in manifest.groupby("gsm", sort=True):
        entries = {r.type: r for r in group.itertuples(index=False)}
        start = time.monotonic()
        count_path = get_file(entries["counts"])
        meta_path = get_file(entries["metadata"])
        result = compare_sums(gsm, count_path, meta_path)
        result["elapsed_seconds"] = round(time.monotonic() - start, 1)
        results.append(result)
        print(gsm, result["status"], result["ec"], flush=True)
        # Never persist or commit downloaded raw cell-level source matrices.
        if os.environ.get("IMACS_QC_KEEP_RAW") != "1":
            count_path.unlink(missing_ok=True)
            meta_path.unlink(missing_ok=True)
    table = []
    for r in results:
        table.append({
            "gsm": r["gsm"], "status": r["status"], "genes": r["genes"],
            "fractional_values_present": r["fractional_values_present"],
            **{"ec_" + k: v for k, v in r["ec"].items()},
            **{"all_" + k: v for k, v in r["all_cells"].items()},
            "ec_nFeature_exact_fraction": r["ec_nFeature_exact_fraction"],
            "ec_nFeature_median_difference": r["ec_nFeature_median_difference"],
        })
    pd.DataFrame(table).to_csv(OUT / "13_GSE225948_full_matrix_normalization_audit.csv", index=False)
    report = {
        "dataset": "GSE225948", "preparations": 12,
        "overall": "PASS" if all(r["status"] == "PASS" for r in results) else "REVIEW",
        "threshold": "At least 95% of annotated EC cells must have full-matrix sum within 5% of nCount_RNA.",
        "interpretation_limit": "No RNA-to-protein, PET, causal, or individual-animal inference.",
        "samples": results,
    }
    (OUT / "13_GSE225948_full_matrix_normalization_audit.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print("OVERALL", report["overall"], "results written to", OUT, flush=True)


if __name__ == "__main__":
    main()
