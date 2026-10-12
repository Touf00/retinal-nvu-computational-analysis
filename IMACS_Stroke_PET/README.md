# IMACS Stroke/NVU PET screen

Temporary isolated working branch for the IMACS Aix-Marseille project screen. Nothing here is merged into the retinal-NVU project on `main`.

## Purpose

Screen cerebral-ischemia molecular-imaging and BBB-delivery candidates without assuming any candidate must win.

Primary order:
1. GSE174574 — replicated 24 h mouse MCAO (3 sham vs 3 MCAO): animal-level pseudobulk.
2. GSE227651 — sham/day 1/day 3/day 7 tMCAO: descriptive temporal trajectories only.
3. GSE233815 — spatial + sc/snRNA + MRI/IHC: spatial adjudication.
4. GSE250245 — rat 48 h MCAO: cross-species/model relevance after reconstructing sample dependence.
5. PET-targetability and BBB-delivery gates only after biology survives.

## /TEAMRED rules

- Cells are not biological replicates.
- One-sample-per-timepoint datasets are never used for inferential p-values.
- RNA abundance is not protein abundance, luminal accessibility, ligand affinity, or PET feasibility.
- Post-stroke uptake can reflect blood flow, BBB leakage, blood volume, metabolism, and clearance—not only target binding.
- Host tracers AP747/APJ and sCD146/AMOT are screened exactly like external candidates.
- RGD/alpha-v-beta-3 and TSPO are comparators, not novelty claims.

## Colab

Open the primary replicated screen:
https://colab.research.google.com/github/Touf00/retinal-nvu-computational-analysis/blob/imacs-stroke-pet-screen/IMACS_Stroke_PET/notebooks/01_primary_GSE174574_screen.ipynb

This branch is a temporary workspace because the connected GitHub action set can create branches/files but cannot create a new repository. Once an empty dedicated repository exists, these files can be moved there unchanged.
