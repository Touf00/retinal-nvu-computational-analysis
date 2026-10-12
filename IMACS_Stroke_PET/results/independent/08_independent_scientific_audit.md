# /TEAMRED validity audit — independent sorted-endothelial and astrocyte-endfoot datasets

- GSE223744: 4 ischemia vs 4 **naïve** animals, CD31+ FACS sorted endothelium. The 24h intervention includes surgery-related changes not controlled by sham surgery.
- Input values `expected_counts` contain fractions. They are not raw integer UMI counts. No DESeq2 or inappropriate integer rounding was performed.
- Exact label permutation n=4/4 means minimum possible two-sided p = 2/70 ≈0.0286, and BH q is across frozen genes only.
- GSE286075: 3 PAIRS contralateral vs ipsilateral endfeet, 6h reperfusion (not 24h).
- GSE286075 uses first (unstranded) STAR gene count column; other strand columns are not independent data.
- Ensembl stable-base ID→symbol map from GSE223744 itself, and conflicts are removed. Unmatched genes indicate mapping limitations, not true absence.
- RiboTag captures ribosome-associated RNA, not protein or a PET signal.
- Rat, whole-brain mouse and sorted EC samples are separate experimental models, not combinable effect-size estimates.
- This is candidate adjudication, not proof of tracer accessibility, specificity, affinity, transport or clinical relevance.
