# Automated /TEAMRED analysis audit

Generated from committed processed results. Execution success is not scientific validation.

| Check | Status | Evidence / limitation |
|---|---|---|
| mouse_acute_samples | PASS | 6 sample rows |
| mouse_acute_animal_level_permutations | PASS | 116 candidate x compartment rows; minimum exact p=0.100; cannot claim p<0.05 |
| mouse_acute_DESeq2 | NOT_COMPLETE | 4 unresolved DESeq2 errors; use exact-label direction analysis only |
| mouse_temporal | PASS | 4 deposited conditions; ONE SAMPLE per condition, descriptive only |
| mouse_temporal_pseudoreplication | PASS | No spurious cell-level inferential p-values |
| rat_48h | PASS | 6 condition/hemisphere matrices; individual rat identities NOT recoverable |
| rat_endothelium_QC | CHECK | Minimum endothelial-labelled nuclei in a group = 12, NOT suitable for vascular validation without improved annotation |
| rat_missing_feature_symbols | CHECK | 4 of frozen genes could not be matched to reference feature symbols; unmatchable is NOT biological absence |
| spatial_sections | PASS | 5 Visium sections: Ctrl x1, day1 x1, day3 x1, day7 x2 |
| spatial_celltypes | CHECK | Mixed Visium spots; no core/periphery annotations; do not assign cell type or lesion location |
| independent_sorted_endothelial_4v4 | PASS | 29 frozen genes, 4 ischemia vs 4 naive control; permutation p min=0.0286; naïve controls not shams |
| independent_astro_endfeet_3paired | PASS | 29 genes with 3 contralateral/ipsilateral pairs; no population-level significance |

## Mandatory interpretation constraints
- Do not pool cells/spots as biological replicates.
- Mouse GSE174574: n=3/group, minimum attainable exact two-sided permutation p=0.1.
- Temporal GSE227651: one sample/condition, trajectories not statistically replicated.
- Rat GSE250245: six group/hemisphere matrices, not independent individual-animal transcriptomes.
- Spatial GSE233814: mixed capture spots and unavailable lesion-region labels in raw count matrices.
- Rat annotation and missing feature symbols require additional review before vascular validation.
- RNA expression cannot establish protein, PET binding, BBB transfer, or human translation.

## Recovery

All outputs are inside `IMACS_Stroke_PET/results/<dataset>/` on the isolated branch.
Each dataset's `source_run_id.txt` identifies the successful GitHub Actions execution.
Raw sequencing inputs are downloaded from GEO and have checksums in per-dataset outputs.
