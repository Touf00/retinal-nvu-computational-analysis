#!/usr/bin/env python3
"""Make the frozen gene × compartment evidence matrix without meta-analysis or pseudoreplication."""
from pathlib import Path
import numpy as np
import pandas as pd

base=Path("IMACS_Stroke_PET/results")
out=base/"integrated"
out.mkdir(parents=True,exist_ok=True)
required=[
base/"primary/06_primary_candidate_exact_permutation_screen.csv",
base/"temporal/02_candidate_timecourse_descriptive.csv",
base/"rat/03_rat_48h_candidate_expression_by_group_and_hemisphere.csv",
base/"spatial/04_spatial_candidate_expression_by_section.csv"
]
missing=[str(x) for x in required if not x.exists()]
if missing:
    (out/"NOT_READY.txt").write_text("Datasets missing; integration blocked, no invented data:\n"+"\n".join(missing))
    print("Integration BLOCKED:",missing)
    raise SystemExit(0)

p=pd.read_csv(required[0]).set_index(["gene","celltype"])
t=pd.read_csv(required[1]).set_index(["gene","broad_celltype","condition"])
r=pd.read_csv(required[2]).set_index(["gene","celltype","prefix"])
z=pd.read_csv(required[3]).set_index(["gene","section"])
# Independent validation cohorts: retain their distinct experimental units and control types.
independent_ec_file=base/"independent/08_independent_FACS_EC_exact_permutation_4v4.csv"
independent_astro_file=base/"independent/08_independent_astro_endfoot_paired_3animals.csv"
ind_ec=pd.read_csv(independent_ec_file).set_index("gene") if independent_ec_file.exists() else None
ind_ast=pd.read_csv(independent_astro_file).set_index("gene") if independent_astro_file.exists() else None

genes=sorted(set(p.index.get_level_values("gene")))
cts=["Endothelial","Astrocyte","Pericyte"]
rows=[]
def fetch(df,idx,col):
    try:
        x=df.loc[idx,col]
        if isinstance(x,pd.Series):return float(x.iloc[0])
        return float(x)
    except (KeyError,TypeError,ValueError):
        return np.nan

def sign(x):
    if pd.isna(x):return "not_available"
    return "higher" if x>0 else ("lower" if x<0 else "unchanged")

for gene in genes:
    for ct in cts:
        effect=fetch(p,(gene,ct),"log2CPM_mean_difference")
        perm=fetch(p,(gene,ct),"exact_label_permutation_two_sided_p")
        d=dict(gene=gene,celltype=ct,
            primary_mouse_log2CPM_difference_24h=effect,
            primary_mouse_exact_two_sided_p=perm,
            primary_mouse_direction=sign(effect),
            primary_mouse_n_animals_per_group=3,
            temporal_independent_animals=0,
            rat_animal_inference="NOT POSSIBLE from deposited pooled group matrices",
            spatial_celltype_inference="NOT POSSIBLE mixed Visium spots")

        sham=fetch(t,(gene,ct,"Sham"),"mean_log1p")
        d["temporal_sham_mean_log1p"]=sham
        for stage in ["MCAO_day1","MCAO_day3","MCAO_day7"]:
            value=fetch(t,(gene,ct,stage),"mean_log1p")
            d[f"temporal_{stage}_minus_sham_mean_log1p"]=value-sham if pd.notna(sham) and pd.notna(value) else np.nan
            d[f"temporal_{stage}_n_cells"]=fetch(t,(gene,ct,stage),"n_cells")
        # Separate source rows do not constitute an independent replication test:
        # temporal dataset contains only one deposited sample per condition.
        d["mouse_24h_vs_day1_reproducibility"]="NOT ESTABLISHED: temporal dataset has no within-timepoint animal replication"

        # Rat severe ipsilateral vs sham right. Recovered compartment counts must pass 50/group
        a=fetch(r,(gene,ct,"sMCAO_ipsi"),"mean_log1p")
        b=fetch(r,(gene,ct,"Sham_R"),"mean_log1p")
        na=fetch(r,(gene,ct,"sMCAO_ipsi"),"n_cells")
        nb=fetch(r,(gene,ct,"Sham_R"),"n_cells")
        d["rat_48h_severe_ipsi_minus_sham_right_mean_log1p"]=a-b if pd.notna(a) and pd.notna(b) else np.nan
        d["rat_48h_qualified_for_descriptive_interpretation"]=bool(pd.notna(na) and pd.notna(nb) and na>=50 and nb>=50)
        d["rat_severe_ipsi_recovered_nuclei"]=na
        d["rat_sham_right_recovered_nuclei"]=nb

        # Tissue spot observations independent of celltype: identical across compartments by definition
        for stage,section in [("control","C1-control"),("day1","B1-D1"),("day3","D1-D3"),
                              ("day7_C1","C1-D7"),("day7_D1","D1-D7")]:
            d[f"spatial_section_{stage}_pct_spots_RNA_detected"]=fetch(z,(gene,section),"pct_spots_detected")
        # Add independent primary biological evidence, not a pooled effect.
        if ct=="Endothelial" and ind_ec is not None:
            d["sorted_CD31_24h_ischemia_vs_naive_log2CPM_difference"]=fetch(ind_ec,gene,"ischemia_minus_naive_log2CPM")
            d["sorted_CD31_4v4_exact_two_sided_p"]=fetch(ind_ec,gene,"permutation_two_sided_p")
            d["sorted_CD31_4v4_BH_q_frozen_panel"]=fetch(ind_ec,gene,"BH_q_across_frozen_gene_panel")
        else:
            d["sorted_CD31_24h_ischemia_vs_naive_log2CPM_difference"]=np.nan
            d["sorted_CD31_4v4_exact_two_sided_p"]=np.nan
            d["sorted_CD31_4v4_BH_q_frozen_panel"]=np.nan

        if ct=="Astrocyte" and ind_ast is not None:
            d["astro_endfoot_6h_paired_ipsi_vs_contra_log2CPM_difference"]=fetch(ind_ast,gene,"mean_paired_ipsilateral_minus_contralateral_log2CPM")
            d["astro_endfoot_6h_n_pairs"]=fetch(ind_ast,gene,"n_pairs")
        else:
            d["astro_endfoot_6h_paired_ipsi_vs_contra_log2CPM_difference"]=np.nan
            d["astro_endfoot_6h_n_pairs"]=np.nan

        d["priority_PET_status"]="NOT RANKED: requires cell-specific protein, accessible epitope, binder and PK"
        rows.append(d)
matrix=pd.DataFrame(rows)
matrix.to_csv(out/"candidate_by_celltype_cross_dataset_evidence.csv",index=False)
summary=matrix[matrix.celltype=="Endothelial"][["gene","primary_mouse_log2CPM_difference_24h",
  "temporal_MCAO_day1_minus_sham_mean_log1p","mouse_24h_vs_day1_reproducibility",
  "rat_48h_qualified_for_descriptive_interpretation","priority_PET_status"]]
summary.to_csv(out/"endothelial_compact_audit.csv",index=False)
(out/"README.md").write_text("""# Frozen cross-dataset candidate synthesis

Each source retains its own quantitative unit. Do not compare effect magnitudes
numerically across independent scRNA, snRNA and Visium experiments.

**Primary (3+3):** exact label-permutation p has theoretical floor 0.1.
**Temporal:** only one deposited sample/condition, descriptive trends only.
**Rat:** six group/hemisphere matrices without recoverable individual rat labels.
**Spatial:** mixed spots, no author region annotations transferred.
**PET feasibility:** not inferred from RNA.

The current cross-study integration deliberately does not calculate 'alignment',
reproducibility percentages, concordance p-values or summary meta-analytic effects.
Sorted CD31 dataset provides an independent n=4/4 ischemia-versus-naive contrast,
not a matched sham-controlled repetition. Endfoot n=3 data are paired ipsilateral
versus contralateral and reflect early ribosome-bound RNA, not global astrocytes.
Time-series data have no within-condition animal replication and marker-based
cell annotations require independent curation before biological validation.
""")
print("Saved",len(matrix),"gene × compartment rows; no statistical meta-analysis")
