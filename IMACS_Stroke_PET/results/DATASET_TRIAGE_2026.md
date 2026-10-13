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

## Human tissue translation — GSE162955 — EXECUTED AND SOURCE-VERIFIED

- Six human postmortem ischemic-stroke donors, **paired infarct core and contralateral tissue**: 12 GPL17586 microarrays.
- Separate notebook `10_human_paired_GSE162955_screen.ipynb` uses original GEO series matrix, authoritative platform annotation and exact **n=6 sign-flip** paired analysis.
- Human **bulk brain tissue** does not supply brain endothelial specificity, BBB accessibility, pericyte specificity, tracer binding or survival-zone signal. Ischemic core cell loss and terminal illness are critical confounders.
- Original downloaded GEO source checksum verified, all six donor/condition labels verified, all probe intensities and 64-permutation paired statistics and BH values independently reproduced. **27/29 panel genes mappable; none reaches frozen-panel human BH q<0.05. Six-donor leave-one-out sign changes in 16/27 genes.** Execution: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38440475042; results: `human_replication_qc/`.
- Any source-file/platform-map failure must be recorded, not hidden by unverified symbol mapping.

## Adult surgical-control human brain FACS endothelium — GSE256490 — EXECUTED

- **Six original TL_patient1–6_EC GEO samples**, GSM8101534–1539; separately downloaded and source-SHA256 fingerprinted. This is an adult temporal lobe control/surgery reference, **not a human AIS case-control study**. Two original source filenames explicitly mention epilepsy.
- After basic per-library QC, **50,251 captured endothelial-enriched cells** from six named patient preparations. All 29 target gene features present; only eight genes detected in >=10% of cells in all six preparations. Cells do not inflate biological sample size.
- Median fraction of cells with detectable RNA across six independent preparations: APLNR **0.12%**, AMOT **4.72%**, MCAM **9.65%**, ICAM1 **16.01%**, PVR **6.12%**, SLC3A2/CD98hc **90.39%**, TFRC **66.88%**, CLDN5 **99.42%**. Low detection does not demonstrate protein absence; high detection does not prove luminal accessibility or BBB transcytosis.
- Executed and archived at https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38441222229; results: `human_vascular_reference/`. Original paper: Wälchli et al. Nature 2024, DOI 10.1038/s41586-024-07493-y.
- Joined as distinct human normal-baseline evidence with acute mouse EC cohorts and six paired human bulk postmortem donors (29-gene translation **hypothesis** gates, not meta-analysis), execution https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38441801735; results: `target_validation_gate/`.

## Other human candidates ruled out as direct AIS cerebral vascular replication

- GSE248760: four sickle-cell patients with versus four without prior ischemic stroke, **peripheral blood-derived endothelial colony-forming cells**. Genuine donor contrast, but not an ischemic human brain microvascular study. Context only.
- GSE331114: poststroke brain immune **CD45-enriched** cells, not unbiased vascular endothelial cells.
- GSE308369: HUVECs stimulated with human stroke plasma, an *in vitro peripheral endothelial model* rather than stroke-patient brain endothelium.

## Exclusions and limits

- **GSE331114** is a public 2026 human AIS brain single-cell study with 4 stroke and 4 nonstroke brain specimens, but the GEO sampling protocol specifically says **CD45+ immune cells were isolated/enriched** from brain. It can inform immune biology; it is **not** an independent unbiased human brain-vascular surface target test.
- **GSE256490 human brain vasculature atlas** is a powerful healthy adult and disease *reference*, but is not automatically a paired human ischemic-stroke dataset. If used for normal off-target/cell-type expression, inspect exact donor diagnoses and keep it separate from ischemia effect inference.
- **SCP3078** 5-GB integrated mouse MCAO atlas expression remains inaccessible from the existing runner route; the portal returned authorization denial. It also reuses existing studies: avoid double-counting its integrated cells as independent validation.
- **GSE250245 rat** pooled six group/hemisphere matrices lack recoverable per-animal counts; reference curation and gene-symbol audit remain necessary before rat vascular conclusions.

## Minimum scientific decision gates

RNA temporal, sample-paired and spatial results are complementary but **must not be pooled into one effect size or spurious meta-analysis**. Prospective PET target selection requires protein presence, in-vivo luminal/accessibility mapping, normal-organ burden, binder/tracer chemistry, nonspecific plasma/BBB-leakage controls, and clinical translation. BBB-shuttle use additionally requires capillary depletion, intact therapeutic cargo verification, parenchymal/cell uptake, knockdown and exposure-to-pharmacology mapping.

**No candidate is 'validated for PET' from these datasets alone.**
