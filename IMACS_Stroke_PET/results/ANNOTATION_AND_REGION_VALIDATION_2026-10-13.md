# IMACS — original cell identity and lesion-region audit, 13 October 2026

Scope: GSE174574, GSE250245, GSE233814. Source-bound limitations; NOT PET target approval.

## GSE174574 — Zheng et al. (PMID 34496660)
Original publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC8721774/ ; GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174574
- Six original GEO libraries (three Sham, three MCAO at 24 hours) re-downloaded with SHA256 matching previously archived input fingerprints.
- Published marker genes Itm2a (endothelium), Hexb (microglia), Aldoc (astrocyte), Plp1 (oligodendrocyte), Ttr (ependymal/CP), Acta2 (VSMC) were present and most highly expressed in the corresponding in-house broad classes on six-sample mean expression.
- Only Itm2a was absent from the original in-house classifier marker panel; the other five markers partially reuse inputs used for classification and cannot be labeled independent validation.
- Ttr detected in essentially all cells across many unrelated in-house broad cell classes, despite highest mean level in Ependymal_CP. Possible ambient/background signal must be audited before interpreting isolated Ttr positivity.
- In the full-order rerun (restoring the original UMAP and raw-layer setup) with Scanpy 1.11.5, 43 sample-by-class count slots differed from archived values, maximum absolute difference 207 cells (MCAO2 Ependymal_CP); original source SHA256 files were identical. The earlier shortened marker run had 47 slots / maximum 52, and is superseded. The changes in broad labels could affect pseudobulk estimates; do not claim exact label or DE reproducibility from this marker audit.
- Original paper resolved 17 types (including arterial, capillary, venous EC). Our pipeline assigns broad computational categories and has NOT recovered original published per-cell identity mappings.
- Executed validation: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38431719568

## GSE250245 — Bormann et al. (PMID 39043661)
Study: https://www.nature.com/articles/s41467-024-50465-z ; GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250245
- Original design: 11 rats (4 Sham, 3 moderate, 4 severe MCAO) with six deposited severity-by-hemisphere matrices. Individual rat-to-barcode membership cannot be reconstructed from our archived inputs.
- Current in-house pipeline retained 73,371 nuclei versus original article's 68,616 across 29 author-defined groups. Our marker-score assignment is NOT an author-label reproduction.
- Four target features absent from released rat gene list in all six matrices: Apln, Amot, Cldn5, Flt1. This cannot be called biological nonexpression.
- Previously 72 missing-feature rows incorrectly passed count-only eligibility. The original rat notebook was repaired and rerun, preserving 1,276 source comparisons: 450 feature-present count-adequate descriptive rows, 650 low-cell-capture rows, and 176 missing-gene-feature rows.
- The in-house endothelial class remains sparse, including only 12 and 14 recovered nuclei in sham hemispheres; no endothelial target validation or rat-level statistical claims.
- Repaired rerun: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38431919460

## GSE233814 / GSE233815 — Zucha et al. (PNAS 2024)
Original dataset: https://data.mendeley.com/datasets/gnb2dsjms2/1 ; GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233814
- Five spatial sections: control, day1, day3, and TWO day7 sections. Sections must NOT automatically be equated with different animals.
- Original DetailedRegionAnnoShort metadata must be recovered from published processed Seurat spatial_seurat_1DP_13_nygen.rds; prior Mendeley API probes received HTTP 403. No original region metadata has yet been accessed or transferred into existing spot expression tables.
- Updated source exporter: scripts/export_zucha_region_annotations.R. New strict importer: scripts/validate_zucha_author_region_crosswalk.py. Both require genuine author fields; mapping day7 by timepoint is refused.
- Synthetic crosswalk self-test passed including rejection of ambiguous day7 mapping: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38432907826 . This is software validation, NOT biological label recovery.
- No original lesion core/periphery evidence, cell-type-specific expression or inferential p-values is established from current spatial archives.

Next genuine source requirements: mouse author barcode-to-label table (if released); rat original author cluster/barcode mapping and animal provenance (if separable); spatial published Seurat RDS and exact five-section mapping. Raw published inputs must not be committed into this repository.

Boundaries: RNA is not protein localization, vascular accessibility, specific tracer binding, BBB transcytosis, or PET viability. All files remain on imacs-stroke-pet-screen; never merge with retinal NVU main.
