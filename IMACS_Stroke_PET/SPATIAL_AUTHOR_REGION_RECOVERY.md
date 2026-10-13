# Run the original-author Zucha spatial region-label transfer

This step is blocked until the published author Seurat RDS is accessible. No lesion region has yet been assigned from source metadata.

Public source (Zucha et al., PNAS 2024): https://data.mendeley.com/datasets/gnb2dsjms2/1

Required author file: spatial_seurat_1DP_13_nygen.rds. Do not substitute GSE233814 raw Visium expression for this metadata, and do not claim expression-derived clusters identify infarct core/periphery.

1. Download the real author RDS through an authorized browser or obtain an equally well-provenanced original barcode+section+DetailedRegionAnnoShort table.
2. In a computer with R and SeuratObject, export all author metadata fields (no counts or image data are exported) using:
    Rscript IMACS_Stroke_PET/scripts/export_zucha_region_annotations.R spatial_seurat_1DP_13_nygen.rds zucha_author_regions.csv
3. Inspect the exported author_sample, condition, author_barcode and candidate author_meta_* columns to identify FIVE original spatial sections, not FOUR timepoint categories.
4. Create a separate manually checked five-row CSV named zucha_section_map.csv with the exact two headers author_sample,gsm, mapping every unique source section label to one of the five GEO IDs GSM7437221 through GSM7437225. Do not invent sample IDs, do not map both day7 samples by Condition alone.
5. After checking identifiers against the original author data, run:
    python IMACS_Stroke_PET/scripts/validate_zucha_author_region_crosswalk.py --author-metadata zucha_author_regions.csv --section-map zucha_section_map.csv
6. Only if strict conditions, barcode joins and sufficient match coverage PASS, inspect the small descriptive results under results/spatial_author_region_qc. Do not treat spots or sections as animal-level replicates.

Safety check verified in GitHub Actions: https://github.com/Touf00/retinal-nvu-computational-analysis/actions/runs/38432907826

The real original RDS has not been acquired by the existing pipeline. The direct Mendeley API returned HTTP 403 to previous GitHub-based source probes. Do not place multi-GB RDS or raw count matrices in the Git repository. If the original file has no unambiguous five-section sample field, obtain author-provided correspondence rather than guessing.

Metadata should only be used to report the author's original DetailedRegionAnnoShort labels. Anatomical core/periphery semantics and PET target access require separate validation.
