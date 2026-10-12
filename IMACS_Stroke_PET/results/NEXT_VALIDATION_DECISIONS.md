# IMACS /TEAMRED — decision checkpoint after dataset processing

Status checked 2026-10-13. **This is not a winning-target declaration.** Read with [RECOVERY_AND_STATUS.md](../RECOVERY_AND_STATUS.md) and [SCIENTIFIC_AUDIT.md](SCIENTIFIC_AUDIT.md). All input checksums, provenance IDs and analysis code are preserved under this branch.

## What has actually been computed

- GSE174574: 3 sham + 3 MCAO mouse cells, compartment-summarized and sample-level pseudobulk candidate counts. Exact balanced animal-label permutations have p-value resolution of **0.10 minimum**. Automated broad lineage labels need independent curation. PyDESeq2 failure documented; **no DESeq2 significance claim**.
- GSE227651: sham/day 1/day 3/day 7 one deposited library each; descriptive candidate trajectories only.
- GSE250245: rat 48-hour data computed, **pooled severity × hemisphere group counts not individual replicates**; low rat EC recovery and unresolved feature names mean **not independent vascular validation**.
- GSE233814: five Visium sections and 29 frozen gene spot values; **mixed spots**, no usable core/periphery labels transferred yet.
- GSE223744: independent FACS CD31+ endothelium, **four ischemia vs four naïve**, exact label permutations; not a sham-surgery contrast.
- GSE286075: astrocyte endfoot RiboTag, three within-animal hemisphere pairs after 2h ischemia + 6h reperfusion; no p<0.05 possible with exact sign flip.
- SCP3078: source metadata investigated, but expression H5AD requires authorization or alternate verified access. **Not analyzed.**

## What the independent endothelial dataset actually argues

Data below are **ischemia minus naïve mean log2CPM** in GSE223744, not protein measurements. Exact permutation p-values and BH correction are **over the frozen 29-gene panel only**. n=4+4, and surgery-vs-naïve confounding persists.

| Target/pathway | GSE223744 RNA change | Exact two-sided p | Interpretation before PET gate |
|---|---:|---:|---|
| Aplnr / APJ | -1.42 | 0.0857 | APJ elevation **not supported** in this independent endothelial experiment |
| Apln / apelin ligand | +4.44 | 0.0286 | Ligand induction does **not** mean receptor induction |
| Amot / angiomotin | -1.58 | 0.0286 | Host AMOT tracer must **not** be preselected |
| Mcam / CD146 | +2.84 | 0.0286 | MCAM is **not** the AMOT target of the soluble-CD146 tracer |
| Ackr1 | +7.43 | 0.0286 | Strong vascular injury biology but erythrocyte/systemic-sink concern |
| Icam1 | +2.21 | 0.0286 | Acute endothelial activation; non-specific systemic inflammation risk |
| Vcam1 | +1.12 | 0.0286 | Endothelial activation; binder/normal-tissue audit still needed |
| Slc22a8 | -2.11 | 0.0286 | BBB-state alteration; not automatically a drug shuttle |
| Tfrc | -1.81 | 0.0286 | Must not claim ischemia-induced BBB shuttle suitability |
| Igf1r | -0.97 | 0.0286 | Same distinction as TFRC |
| Slc3a2/CD98hc | +0.03 | 0.8857 | Relatively preserved RNA, not proof of productive CNS delivery |

GSE174574 endothelial exploratory pseudobulk direction is also negative for **Aplnr** and **Amot** in the frozen panel. This is more informative than searching only for positive host-tracer evidence, but the two datasets are different models/controls with different measurement technologies; **do not pool effect sizes**.

GSE286075 detects lower **Maob** ribosome-associated RNA in 3/3 ischemic astrocyte-endfoot samples at **6h reperfusion**, whereas the GSE227651 *whole astrocyte* descriptive readout suggested potentially higher Maob RNA at **day 3**. These are not contradictory tests of one common estimand: the times, tissue compartments and technologies differ. Neither result proves MAO-B protein nor astrocyte PET validity.

## Targets not yet ready for radiotracer or shuttle selection

No candidate earns a PET/BBB-delivery GO from RNA datasets alone. The next essential filters are:
1. **Human and rat brain injury protein-level expression** and localization.
2. Luminal accessibility for intravascular vascular-PET targets; neuronal/abluminal penetration for parenchymal tracers.
3. Binding ligand, competition/blocking, tracer plasma metabolism and isotope suitability.
4. Normal organ sinks and expected specific signal relative to blood pool/perfusion/leakage.
5. Direct **intact therapeutic**, capillary-depletion and target-engagement assays for BBB delivery proposals, not only labeled-brain uptake.

## Recovery and next tasks

Priority A — Obtain the **author's detailed spatial region annotation** for GSE233814 (published processed Seurat Mendeley: DOI 10.17632/gnb2dsjms2.1); current results cannot support lesion core/penumbra localization. Original paper: Zucha et al., PNAS 2024.

Priority B — Recover original annotated cell identities for acute mouse and rat, then rerun target expression under reference-validated labels. The very low rat endothelial count makes current cell-restricted comparisons unreliable.

Priority C — Determine whether SCP3078 is accessible through an authorized download path. Do **not** label the atlas complete or claim independent replication merely by concatenating overlapping source studies.

Priority D — Run a separate **candidate tracer and BBB shuttle targetability** screen against human proteins, binding chemistry, normal organ expression and drug-distribution controls.

Every result file under this branch is intentionally kept separate from the original retinal `main` branch; do not merge branches. No raw GEO files or confidential industry partner data are archived.
