# Frozen cross-dataset candidate synthesis

Each source retains its own quantitative unit. Do not compare effect magnitudes
numerically across independent scRNA, snRNA and Visium experiments.

**Primary (3+3):** exact label-permutation p has theoretical floor 0.1.
**Temporal:** only one deposited sample/condition, descriptive trends only.
**Rat:** six group/hemisphere matrices without recoverable individual rat labels.
**Spatial:** mixed spots, no author region annotations transferred.
**PET feasibility:** not inferred from RNA.

The current cross-study integration deliberately does not calculate 'alignment',
reproducibility percentages, concordance p-values or summary meta-analytic effects.
Time-series data have no within-condition animal replication and marker-based
cell annotations require independent curation before biological validation.
