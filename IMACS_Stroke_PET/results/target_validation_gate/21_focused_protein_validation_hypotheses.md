# IMACS focused protein-validation hypotheses — NOT PET nominations

Frozen panel: all 29 genes remain in the machine-readable evidence table. The following nine are hypothesis-focused laboratory comparators (six vascular imaging hypotheses, two delivery handles, one astrocyte target), not genome-wide statistically prioritized winners.

Source: six GSE174574 mouse 24h MCAO/Sham libraries; eight GSE223744 FACS sorted mouse ischemia/naive libraries; GSE225948 pooled stroke/sham 4+4 with batch concerns; six GSE162955 paired HUMAN postmortem whole-tissue donors; SIX independently labeled human adult temporal-lobe control FACS brain EC samples from GSE256490. Different estimands MUST NOT be meta-analyzed.

**Human six-paired infarct-core whole-tissue dataset has ZERO genes with q<0.05. Mouse exact 3v3 permutation p-values have floor 0.10. No candidate is PET-ready.**

| Gene | Track | Source-aware caution | Next discriminating experiment |
|---|---|---|---|
| MCAM | surface-vessel/mural localization hypothesis | CD146 different from sCD146 ligand AMOT; human EC baseline variable | capillary/venule vs pericyte/VSMC protein, normal-organs, ligand cross-reactivity |
| ICAM1 | vascular inflammation surface hypothesis | broad organ endothelial activation, human core statistics not significant | human brain ischemia EC surface IHC and specific blood-pool/BBB-corrected binder |
| PVR | immune-endothelial interface hypothesis | human baseline EC low to intermediate; no validated brain PET binder here | human endothelial surface protein and nonvascular/background binding |
| VCAM1 | vascular inflammation comparator hypothesis | human postmortem bulk direction differs from mouse acute; normal systemic inflammation | species-specific brain EC surface protein, timepoint and organ sink |
| APLNR | host APJ falsification hypothesis | mouse acute and sorted EC direction negative, human normal EC RNA sparse | luminal APJ protein after ischemia, receptor blocking vs ischemic blood delivery |
| AMOT | host sCD146 binding-target falsification hypothesis | human normal EC RNA low; AMOT surface fraction and accessibility unresolved | human brain ischemic endothelial accessible AMOT epitope, sCD146 competition and autoradiography |
| SLC3A2 | BBB delivery platform comparator | strong human normal EC RNA; transporter expression is NOT intact cargo transcytosis | construct-specific exposure, capillary depletion and target engagement |
| TFRC | BBB delivery platform comparator | substantial normal human EC expression but mouse ischemic RNA decreases | TfR shuttle construct binding, capillary depletion and drug effect |
| MAOB | astrocytic PET biology comparator | human EC expression limited; astrocyte modality, not vascular-surface marker | reactive astrocyte-specific enzyme/protein activity and tracer specificity |

## Evidence flags for the nominated protein-hypothesis set

| Gene | Human normal EC median RNA detection | Mouse EC 2-cohort direction | Human whole tissue 6-pair mean direction | Human donor leave-one-out sign changes |
|---|---:|---|---|---:|
| MCAM | 9.65% | Concordant | increase | 0 |
| ICAM1 | 16.01% | Concordant | increase | 1 |
| PVR | 6.12% | Concordant | increase | 0 |
| VCAM1 | 11.28% | Concordant | decrease | 1 |
| APLNR | 0.12% | Concordant | decrease | 0 |
| AMOT | 4.72% | Concordant | increase | 1 |
| SLC3A2 | 90.39% | Discordant | increase | 1 |
| TFRC | 66.88% | Concordant | decrease | 0 |
| MAOB | 3.88% | Concordant | decrease | 0 |

## Scientific decision gate

- Low transcript detection in control human EC may reflect true low baseline, sampling/capture dropout, anatomical environment or regulation. It cannot prove protein absent or rule out induction after stroke.
- FACS-sorted human controls are surgical temporal-lobe samples, including epilepsy-referenced filenames; not baseline for all brain regions nor stroke controls.
- Human stroke six-donor core vs contralateral is whole-brain postmortem tissue, affected by cell death and heterogeneous composition. No endothelial-specific paired ischemic effect is established.
- GSE225948 EC2 fraction is highly shifted in pooled stroke vs sham preparations with unexchangeable batch dates; not independent disease attribution.
- Protein localization, accessible epitope, blocking, competition, radiochemistry, physiological transport, normal-organ sink, BBB leak and perfusion confounding are mandatory separate stages.
- **Only after a small protein-screen shortlist survives source- and vessel-validated ischemic localization should the original Zucha spatial DetailedRegionAnnoShort metadata be used to test region specificity.**

Method: no statistical meta-analysis, no p-value combination, no whole-tissue-to-cell-type causal claims. Human tissue direction matching mouse EC is exploratory only.
