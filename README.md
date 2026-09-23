# Retinal endothelial Norrin–BRB computational analysis

Reproducible secondary analysis of three public mouse retinal single-cell RNA-sequencing resources. The repository tests whether a frozen Norrin/β-catenin-associated endothelial module covaries with a frozen blood–retinal barrier (BRB) transcriptional module across complementary datasets.

This is an association analysis, not a pathway-activity assay, causal test, clinical-biomarker study, or genotype-level replication analysis.

## Headline results

| Dataset | Analysis unit | Clean ECs | Validated result |
|---|---:|---:|---|
| MRCA | 11 sufficiently sampled sequencing cohorts | 1,436 | Descriptive random-effects `r = 0.206` (95% CI `0.099–0.309`), `I² = 66.5%` |
| Furtado, GSE282775 | 5 pooled genotype libraries | 3,017 | Positive within-library Norrin–BRB Spearman correlations (`ρ = 0.122–0.236`), essentially unchanged after AV adjustment |
| Zarkada, GSE175895 | 3 WT + 3 Alk5 libraries | 683 WT + 527 Alk5 | D-tip-like > S-tip-like for Norrin, BRB, junction, and transport modules in 6/6 libraries |

The MRCA heterogeneity is substantial and should be reported, not averaged away. Furtado libraries are pooled preparations, and Zarkada state-specific cell counts are unequal. Individual cells are not independent biological replicates; all within-library cell-level comparisons are descriptive.

## Frozen modules

- Primary Norrin-associated module: `Fzd4`, `Lrp5`, `Tspan12`, `Lef1`
- BRB module: `Cldn5`, `Ocln`, `Tjp1`, `Mfsd2a`, `Slc2a1`, `Abcb1a`
- Junction sensitivity module: `Cldn5`, `Ocln`, `Tjp1`
- Transport sensitivity module: `Mfsd2a`, `Slc2a1`, `Abcb1a`
- Legacy Norrin sensitivity module only: `Ndp`, `Fzd4`, `Lrp5`, `Tspan12`, `Ctnnb1`

The legacy definition never replaces the primary module.

## Repository contents

- `notebooks/retinal_nvu_computational_analysis.ipynb` — linear download, QC, analysis, validation, and figure pipeline
- `results/reference/` — compact frozen outputs from the validated execution
- `results/reference/figure_computational_triangulation.png` — validated three-panel reference figure
- `results/manifest.json` — machine-readable data, module, result, caveat, and provenance manifest
- `requirements.txt` and `environment.yml` — environment specifications
- `tests/test_release_package.py` — dependency-free structural and reference-result checks
- `ZENODO_RELEASE.md` — exact GitHub-to-Zenodo release procedure

Large public inputs and generated cell-level tables are deliberately excluded.

## Public datasets

| Resource | Record | Input |
|---|---|---|
| Mouse Retina Cell Atlas (MRCA) | [Zenodo 10815031](https://doi.org/10.5281/zenodo.10815031) | 3.6 GB all-cells `.h5ad`, MD5 `cbbd7cc9802330a7002946757b995f26` |
| Furtado et al. | [GEO GSE282775](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282775) | Five processed 10x matrices in `GSE282775_RAW.tar` |
| Zarkada et al. | [GEO GSE175895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175895) | Three WT and three Alk5 processed 10x matrices in `GSE175895_RAW.tar` |

The source datasets remain governed by their own repository records and terms; the MIT license in this repository applies to this analysis code, not to redistributed source data. No raw source data are redistributed here.

## Run from a clean environment

Python 3.13 was used for the validated run. Core observed package versions are pinned; packages whose exact versions were not recorded in the original execution are bounded in `requirements.txt`. The notebook contains fail-fast assertions for cell counts, cluster structure, and headline results to expose version drift.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python -m ipykernel install --user --name retinal-nvu --display-name "Python (retinal-nvu)"
```

To let the notebook download the public inputs, set `RETINAL_NVU_DOWNLOAD=1` before launching Jupyter:

```bash
# Windows PowerShell
$env:RETINAL_NVU_DOWNLOAD = "1"
jupyter lab

# macOS/Linux
export RETINAL_NVU_DOWNLOAD=1
jupyter lab
```

Alternatively, place files under:

```text
data/raw/
├── MRCA/MRCA_all_cells.h5ad
├── GSE282775/<five extracted 10x triplets>
└── GSE175895/<six extracted 10x triplets>
```

Run every notebook cell from top to bottom. Generated outputs go to `results/generated/`, which is ignored by Git.

Expected resources: approximately 8–12 GB of download/storage headroom and at least 16 GB RAM. The Zarkada matrices contain millions of barcodes and dominate runtime and memory pressure.

## Reproducibility status

The compact reference files were transcribed from a successful executed notebook whose SHA-256 is recorded in `results/manifest.json`. During repository hardening, the final pipeline was statically compiled cell-by-cell and the release tests were run. A new full raw-data rerun was not performed during packaging because the multi-gigabyte public inputs are intentionally not bundled. Therefore the honest release gate is: **run the cleaned notebook once in the declared environment and require every assertion to pass before tagging the archival release**.

## Interpretation boundary

Use language such as “Norrin-associated endothelial identity covaries with BRB-associated transcriptional identity.” Do not claim that these module scores directly measure Norrin pathway activity, prove that Norrin causes the BRB state, establish genotype-level effects, or validate a clinical biomarker. Shared developmental or maturation programs remain plausible explanations.

## Citation and authorship metadata

`CITATION.cff` currently uses an entity-level contributor label because personal author names and ORCIDs were not supplied during packaging. Replace that entry with the actual software authors before the first public release. Do not mint a DOI with placeholder authorship.

## License

Analysis code is offered under the MIT License. Confirm that this matches all contributors’ intentions before release.

