import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_csv(name):
    with (ROOT / "results" / "reference" / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


required = [
    "README.md",
    "LICENSE",
    "CITATION.cff",
    "requirements.txt",
    "environment.yml",
    ".gitignore",
    "ZENODO_RELEASE.md",
    "notebooks/retinal_nvu_computational_analysis.ipynb",
    "results/manifest.json",
]
for relative in required:
    assert (ROOT / relative).is_file(), f"Missing release file: {relative}"

notebook_path = ROOT / "notebooks" / "retinal_nvu_computational_analysis.ipynb"
notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
assert notebook["nbformat"] == 4
assert len(notebook["cells"]) >= 15
source = "\n".join(str(cell.get("source", "")) for cell in notebook["cells"])
for forbidden in ["/content/", "RECOVERED_TARGETS", "capillary_vs_venule_z_example", "P6_tip_barrier_mean_example"]:
    assert forbidden not in source, f"Obsolete/private token remains: {forbidden}"
for required_gene in ["Fzd4", "Lrp5", "Tspan12", "Lef1", "Cldn5", "Ocln", "Tjp1", "Mfsd2a", "Slc2a1", "Abcb1a"]:
    assert required_gene in source

mrca = json.loads((ROOT / "results" / "reference" / "MRCA_random_effects_summary.json").read_text(encoding="utf-8"))["PRIMARY_NORRIN_vs_BRB"]
assert mrca["k"] == 11
assert math.isclose(mrca["pooled_r"], 0.206009731473561, abs_tol=1e-12)
assert math.isclose(mrca["I2_percent"], 66.54151894364368, abs_tol=1e-10)

furtado = read_csv("Furtado_5library_Norrin_BRB_summary.csv")
assert len(furtado) == 5
assert sum(int(row["n_cells"]) for row in furtado) == 3017
assert all(float(row["raw_Norrin_vs_BRB_rho"]) > 0 for row in furtado)
assert all(abs(float(row["change_after_AV_adjustment"])) <= 0.004 for row in furtado)

zarkada = read_csv("Zarkada_FINAL_6library_S_tip_D_tip_summary.csv")
assert len(zarkada) == 6
for module in ["norrin", "barrier", "junction", "transport"]:
    column = f"delta_D_minus_S_{module}_score"
    assert all(float(row[column]) > 0 for row in zarkada)

manifest = json.loads((ROOT / "results" / "manifest.json").read_text(encoding="utf-8"))
assert manifest["reference_results"]["mrca"]["pooled_r"] == mrca["pooled_r"]
assert manifest["reference_results"]["furtado"]["clean_ecs"] == 3017
assert manifest["reference_results"]["zarkada"]["directionally_positive_libraries"] == 6

print("Release package checks passed.")
