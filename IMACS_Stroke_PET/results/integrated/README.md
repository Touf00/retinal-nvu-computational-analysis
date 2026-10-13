# Frozen cross-dataset candidate synthesis

Each source retains its own quantitative unit. Do not compare effect magnitudes
numerically across independent scRNA, snRNA and Visium experiments.

**Primary (3+3):** exact label-permutation p has theoretical floor 0.1.
**Temporal:** only one deposited sample/condition, descriptive trends only.
**Rat:** six group/hemisphere matrices without recoverable individual rat labels. A candidate needs BOTH a present deposited feature and >=50 cells in each contrasted compartment even for descriptive inclusion; no animal-level differential expression can be inferred.
**Spatial:** mixed spots, no author region annotations transferred.
**PET feasibility:** not inferred from RNA.

The current cross-study integration deliberately does not calculate 'alignment',
reproducibility percentages, concordance p-values or summary meta-analytic effects.
Sorted CD31 dataset provides an independent n=4/4 ischemia-versus-naive contrast,
not a matched sham-controlled repetition. Endfoot n=3 data are paired ipsilateral
versus contralateral and reflect early ribosome-bound RNA, not global astrocytes.
Time-series data have no within-condition animal replication and marker-based
cell annotations require independent curation before biological validation.
