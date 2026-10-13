# Dataset triage after the 2026-10-13 execution round

## Independent 2026 spatial tMCAO — GSE337578 — EXECUTED

- Public mouse 60-min tMCAO, day 3, **one GEO spatial library** (GSM9858163).
- 30,560 spatial barcodes, 32,286 feature records; candidate RNA values and matching coordinates processed successfully in GitHub Actions run **38325519991**.
- GEO spatial columns **Cell_Barcode, X, Y** are audited and matched. Spot-level expression maps are generated.
- 29 frozen genes queried; candidate detection rates are **spatial-barcode frequencies**, *not* cell-type expression or comparable proportions to scRNA-seq.
- No biological replicates, no sham contrast, no validated lesion core/peri-infarct annotation. No statistical DE, differential cell abundance, or 'core enriched' claim.
- Outputs archived under `results/spatial_day3/` with checksum and source run ID. Figures are preserved in the successful Actions artifact.

## Original Zucha spatial regions — GSE233814 / GSE233815 — PARTIALLY BLOCKED

- Five Visium sections and 29 frozen-gene spot counts are computed and archived.
- The author's `DetailedRegionAnnoShort` is encoded in the published processed Seurat RDS (Mendeley DOI: **10.17632/gnb2dsjms2.1**). It is not reconstructed reliably from unannotated GEO count matrices.
- Public Mendeley API variants all returned HTTP **403 Cloudflare** to the GitHub runner (verified in the source-probe archive). Do not bypass access controls or infer labels from gene expression.
- The associated processed Seurat objects include `seurat_1stSpatial.rds` (~148 MB), `seurat_2ndSpatial.rds` (~100 MB), and `spatial_seurat_1DP_13_nygen.rds` (~3.01 GB). The **3.01 GB final object** is the verified potential source of original detailed lesion-region labels.
- The paper's open interactive Nygen ScarfWeb may also display region annotations. It is a visual aid, not a proven machine-readable export.
- When an authorized user download of the processed object or an author-exported **barcode + sample ID + DetailedRegionAnnoShort CSV** becomes available, join on both specimen and barcode; do not use barcode alone across sections. Hold out day-7 section replication explicitly.

## Human tissue translation — GSE162955 — INITIATED

- Six human postmortem ischemic-stroke donors, **paired infarct core and contralateral tissue**: 12 GPL17586 microarrays.
- Separate notebook `10_human_paired_GSE162955_screen.ipynb` uses original GEO series matrix, authoritative platform annotation and exact **n=6 sign-flip** paired analysis.
- Human **bulk brain tissue** does not supply brain endothelial specificity, BBB accessibility, pericyte specificity, tracer binding or survival-zone signal. Ischemic core cell loss and terminal illness are critical confounders.
- Any source-file/platform-map failure must be recorded, not hidden by unverified symbol mapping.

## Exclusions and limits

- **GSE331114** is a public 2026 human AIS brain single-cell study with 4 stroke and 4 nonstroke brain specimens, but the GEO sampling protocol specifically says **CD45+ immune cells were isolated/enriched** from brain. It can inform immune biology; it is **not** an independent unbiased human brain-vascular surface target test.
- **GSE256490 human brain vasculature atlas** is a powerful healthy adult and disease *reference*, but is not automatically a paired human ischemic-stroke dataset. If used for normal off-target/cell-type expression, inspect exact donor diagnoses and keep it separate from ischemia effect inference.
- **SCP3078** 5-GB integrated mouse MCAO atlas expression remains inaccessible from the existing runner route; the portal returned authorization denial. It also reuses existing studies: avoid double-counting its integrated cells as independent validation.
- **GSE250245 rat** pooled six group/hemisphere matrices lack recoverable per-animal counts; reference curation and gene-symbol audit remain necessary before rat vascular conclusions.

## Minimum scientific decision gates

RNA temporal, sample-paired and spatial results are complementary but **must not be pooled into one effect size or spurious meta-analysis**. Prospective PET target selection requires protein presence, in-vivo luminal/accessibility mapping, normal-organ burden, binder/tracer chemistry, nonspecific plasma/BBB-leakage controls, and clinical translation. BBB-shuttle use additionally requires capillary depletion, intact therapeutic cargo verification, parenchymal/cell uptake, knockdown and exposure-to-pharmacology mapping.

**No candidate is 'validated for PET' from these datasets alone.**
