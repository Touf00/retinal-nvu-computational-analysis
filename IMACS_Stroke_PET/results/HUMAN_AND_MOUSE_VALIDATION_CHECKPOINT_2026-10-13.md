# IMACS human cerebral vascular evidence and replicate-triage checkpoint (13 Oct 2026)

All analyses and outputs below were run on `imacs-stroke-pet-screen` without touching `main`. No candidate is validated for PET or BBB transcytosis by transcriptomics.

## Human ischemic brain tissue — GSE162955
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162955
- Six deceased ischemic stroke patients, matched infarct-core vs contralateral brain homogenates (12 GPL17586 processed arrays). No separation of brain EC from parenchymal cells or vasculature.
- Independent original GEO source and six pair titles/source SHA256 checked, original 27-platform-mappable gene probe intensities reconstructed, exact paired 64-sign-flip tests and frozen-panel BH recovered identically.
- APELA and ACKR1 lacked unambiguous original microarray platform probe mappings. This does not mean they are absent biologically.
- No human gene reached adjusted BH q<0.05 among 27 testable frozen genes; 16/27 reverse mean sign in some one-patient-omitted sensitivity. Core tissue necrosis/cellular composition/terminal illness can dominate bulk expression.
- Execution: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38440475042
- Outputs: results/human_replication_qc/ (executed and archived).

## Human source-specific baseline — GSE256490, Wälchli et al. Nature 2024
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256490
- Authors: https://pmc.ncbi.nlm.nih.gov/articles/PMC11324530/ ; human vascular atlas includes other brain diseases but NOT an AIS brain disease contrast.
- Six specific adult TL control patient-named FACS endothelial source libraries GSM8101534 to GSM8101539, all backed by original source metadata and per-sample GEO matrices. Two matrix filenames mention Epilepsy; these are temporal lobe surgical specimens, NOT a healthy population and NOT stroke.
- Six original GEO matrices independently re-downloaded, 50,251 source cells retained after a documented basic QC pass; original SHA256, sample-level and marker results archived.
- All 29 frozen candidate symbols were present as source features in all six samples; this alone does not imply transcript detection. 8 of 29 genes had detection in at least 10% of captured cells in all six source preparations.
- Median per-preparation detection fraction: APLNR 0.12%, AMOT 4.72%, MCAM 9.65%, ICAM1 16.01%, PVR 6.12%, SLC3A2 90.39%, TFRC 66.88%, IGF1R 43.96%, CLDN5 99.42%, MFSD2A 80.61%.
- Per-cell normalized RNA baseline is descriptive only, is detection-depth- and QC-sensitive, and cannot prove protein in vessels, blood-facing epitope, delivery or PET.
- Execution: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38441222229 ; outputs: results/human_vascular_reference/.

## Mouse EC biological replication (distinct disease models and controls)
- GSE174574: three sham vs three MCAO mouse libraries, 24 h, marker-scored broad EC class. Exact balanced p-value floor 0.10. All 29 EC gene effect signs retained under in-house reannotation sensitivity, without confirming original author cell ID classes.
- GSE223744: four ischemic mice vs four NAIVE controls, sorted CD31+ brain EC, 24 h. Exact 4v4 p-value floor 0.0286; q values over frozen panel only. Naive-v-surgery group confounding blocks causal clean stroke attribution.
- Of 29 frozen candidate genes 25 have the same effect direction across these two mouse endothelial cohorts; this is descriptive direction comparison, not quantitative effect meta-analysis or independent human replication.
- GSE225948 day2 sham-v-stroke FACS preparation comparisons have 4 pooled preparations per arm, severe captured EC2 compositional shift and processing/batch non-exchangeability; no 28-panel BH significance. Do not count as equivalent to 8 independent animals.

## New or alternative human discovery leads (scientific inclusion gates)
- GSE248760 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248760: n4 sickle-cell anemia with prior ischemic stroke versus n4 without stroke, **blood-derived endothelial colony-forming cells**. Genuine human donor replication, but peripheral progenitor identity, sickle-cell disease and past event confound acute CNS vascular targeting; context only, not independent normal adult brain EC AIS replication.
- GSE331114 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331114: human AIS n4 and nonstroke n4 brain parenchyma with CD45 immune enrichment; useful immune response, NOT unbiased human endothelial transcriptome.
- GSE308369 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308369: human umbilical-vein endothelial culture exposed to plasma from stroke patients; stimulation model, NOT cerebral microvascular specimens or paired stroke-brain tissue.
- GSE256490 is valuable to model normal human cerebrovascular transcriptional background; tumor and AVM subcohorts should not be mislabeled stroke.

## Integrated 29-gene gate and priority experimental tracks
- Run: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38441801735
- Full table: results/target_validation_gate/21_frozen_29_gene_human_mouse_evidence_gates.csv
- Protein-hypothesis tracks: MCAM, ICAM1, PVR, VCAM1 as vascular surface context; APLNR and AMOT as host-tracer falsification controls; SLC3A2 and TFRC as BBB-cargo shuttle comparators; MAOB separate astrocyte-associated PET biology. This is NOT PET-ready target ranking.
- Twelve of 27 measurable human postmortem whole-tissue gene effects shared direction with both mouse EC cohorts. Zero reached human q<0.05. Do not combine effect sizes or infer human vascular induction.

## Stop/go
Mandatory next external validation is human ischemic brain surface protein and lesion-region localization versus normal-organ endothelial/pericyte signal, binder competition, tracer metabolism and blood-pool/perfusion/BBB leakage corrections. For cargo delivery demonstrate intact parent/cargo and parenchymal uptake via capillary depletion and pharmacodynamic target engagement.
Original Zucha `DetailedRegionAnnoShort` spot metadata is still awaited from author or processed RDS; do not guess labels. A focused protein screen is more important than processing a 3 GB RDS solely for authenticity.
