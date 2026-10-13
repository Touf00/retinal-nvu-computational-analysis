# Independent sham-controlled stroke endothelial screen: GSE225948

## Why this dataset was added

The first IMACS screening stack emphasized acute single-cell sampling and temporally unreplicated atlases. GSE225948 (Garcia-Bonilla et al., PMID 38177281) provides a **separate** brain-cell experiment with four young-male operated-sham preparations and four young-male tMCAO preparations **both at day 2**. This is a suitable orthogonal candidate challenge, not automatically a stroke vessel-PET validation.

## Original GEO design and actual unit of replication

- GEO samples GSM7060815–GSM7060818: day-2 sham, four deposited preparations.
- GSM7060819–GSM7060822: day-2 stroke, four deposited preparations.
- All deposited samples are male C57BL/6J young brain samples in this planned analysis.
- GEO protocols describe **four to five pooled hemispheres for each experiment**; the unit is a pooled preparation. Do not call this four independent animals per group.
- Samples are FACS-selected for brain immune/myeloid and endothelial populations, not an unbiased whole-brain cell mixture.
- Day-14 stroke samples **lack matched day-14 sham**; do not infer d14 injury effects from simple d2 stroke comparison.
- Report subtype-mix changes independently from per-cell transcript change.

## Source formatting audit completed

GitHub Actions source-probe run 38328531278 inspected all eight exact GEO sample filenames and author cell metadata. Count-matrix headers and barcodes match sample-specific metadata layouts, with **genes in rows, cell barcodes in columns**.

The deposited `*_counts.csv.gz` matrices contain *fractional-valued entries*, so analysis must not presume genuine integer UMI counts or feed them directly to negative-binomial DESeq2. Exact balanced-label randomization is run on one normalized value per independent pooled preparation.

Author `metadata.csv` has `parent` and `sub.celltype`; the `parent == EC` identity should be used before any unsupervised candidate mining. Observed EC cells by sample in GSM order:

| Sham | Author-annotated ECs |
|---|---:|
| GSM7060815 | 546 |
| GSM7060816 | 617 |
| GSM7060817 | 321 |
| GSM7060818 | 295 |

| Stroke | Author-annotated ECs |
|---|---:|
| GSM7060819 | 360 |
| GSM7060820 | 458 |
| GSM7060821 | 502 |
| GSM7060822 | 1171 |

**Source annotation shift:** sham samples are dominated by author `EC1` cells (259, 358, 216, 200), while all stroke preparations show much more author `EC2` (170, 271, 193, 383) than sham samples (30, 18, 5, 1). This is a change in *captured annotated-cell composition*, not proof of tissue frequency or transdifferentiation.

## Analysis execution

Notebook: `notebooks/11_GSE225948_FACS_EC_4v4_day2.ipynb`.

Workflow: `.github/workflows/imacs_GSE225948_day2.yml`.

Status at writing: workflow launched; results not yet scientifically audited.

## Required /TEAMRED checks after execution

1. Ensure eight exact sample/barcode joins and author EC parent labels; no inferred cell identities.
2. Acknowledge 4 biological **pooled preparations**, not 4 independent mice, per group.
3. Audit whether stroke effects result from EC subtype composition rather than within-state transcription; compare within EC1 where n >= 10 for all eight samples. EC2 Sham n1 in one prep makes all-eight within-EC2 inference invalid.
4. Report source count matrix's fractional values and normalization conventions; effect is a difference in sample-mean log1p normalized RNA, not transcriptional fold-change or protein.
5. Use exact group-label permutation over 70 assignments; minimum two-sided p = 2/70. Control multiplicity only within the frozen panel, and state that limitation.
6. PET requires accessible protein, ligand binding, BBB/perfusion controls and radiochemical validation; BBB delivery requires capillary depletion/intact drug/target engagement.
