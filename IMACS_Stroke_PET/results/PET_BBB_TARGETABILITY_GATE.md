# Molecular PET and BBB delivery candidate adjudication — /TEAMRED gate (13 Oct 2026)

**Not a target ranking or tracer nomination.** Based on mouse/rat/human RNA analyses already archived in this isolated branch plus cited tracer/BBB-delivery papers. Each row identifies an experiment that would change the decision. Disease relevance, assay feasibility, translational pharmacology and industrial fit must be weighed **separately** from existing host ownership.

## Stroke vascular PET

| Candidate | Existing evidence | RNA challenge from our source-controlled analyses | PET/development red flags | Required stop/go experiment |
|---|---|---|---|---|
| **APJ / APLNR** — host Ga-AP747 | Original host ligand/ischemia PET exists, cerebral ischemia apelin biology plausible | GSE223744 independent FACS EC ischemia-v-naïve Aplnr **-1.42 mean log2CPM**, exact p=.0857; new GSE225948 sham-controlled day2 EC prep pooled mean change +.0137 in mean log1p normalized RNA, p=1.00; GSE174574 acute provisional EC direction also not positive | BBB leakage, perfusion, cell-type mixture and sparse transcript detection; peripheral vascular tracer success does not transfer automatically to cerebral vasculature | In rat/mouse stroke, quantify **APJ protein on blood-facing cerebrovascular endothelium** with controls and autoradiography/competition; PET specificity vs flow/BBB-leakage; compare animal-day kinetics |
| **AMOT / angiomotin** — host Ga-sCD146 | Host PET in hindlimb and myocardial ischemia; **sCD146 tracer binds AMOT**, not MCAM/CD146 | GSE223744 EC Amot **-1.58**, exact p=.0286; GSE225948 day2 EC -0.0539, p=.0571 (panel FDR q=.178); GSE174574 acute EC negative exploratory direction | AMOT topology, accessible ligand-binding fraction and normal-organ distribution unclear; expression changes may reflect vascular subtype composition | Validate extracellular/luminal AMOT epitope under intact BBB and after stroke, blocking with sCD146, timepoint protein and blood-pool controls |
| **MCAM/CD146** | Biologically relevant vascular marker and receptor, but **NOT** the AMOT ligand target detected by Ga-sCD146 PET | GSE223744 MCAM +2.84; GSE225948 author-annotated day2 EC +0.076, p=.0286; FDR q=.133 | Pericyte/VSMC and extracerebral endothelial expression; no assumption of EC specificity; binder and label PK must be checked | Histology on luminal EC vs mural cells + normal-organ uptake/surface proteome before PET chemistry |
| **ICAM1 / VCAM1** | Accessible inflammatory EC response molecules with existing targeting chemistry in other models | Independent EC study ICAM1/VCAM1 increased; day2 GSE225948 ICAM1 +0.071 (p=.0286; q=.133), VCAM1 no clear contrast (p=.886) | Systemic inflammation and normal-organ endothelial signal; target induction can reflect acute injury, not reparative angiogenesis | Brain EC surface protein and perfusion/BBB-corrected tracer specificity; compare target/background with conventional vascular inflammation imaging |
| **ACKR1** | Acute inflammatory venular biology, possible ligand accessibility | Independent ischemic endothelial RNA strongly increased; GSE225948 pooled EC day2 +0.151, p=.0286, q=.133 | Duffy antigen on RBCs creates potential blood pool/erythrocyte binding sink; binder competition and species differences | Demonstrate endothelial-selective accessible ligand without RBC binding, ex-vivo plasma/whole-blood PK |
| **MAO-B** astrocyte-associated PET | Existing human stroke astrocytosis PET with THK5351; modern reversible MAO-B tracers exist | Endfoot ribosome-associated Maob lower at acute 6h in n=3 paired ischemia; one-timepoint-per-condition mouse temporal later day3 appears higher; neither measures protein | Expression/time/compartment mismatch; MAO-B not perfectly selective for all astrocyte states; availability of modern fluorine-18 tracer unknown | Protein/activity mapping across 6h, 1d, 3d, 7d, 14d; compare independent astrocyte reactive state IF with PET kinetics |

## BBB cargo delivery (different screening universe)

| Handle | In-vivo platform maturity | Transcriptome signal limitation | Drug-development decision |
|---|---|---|---|
| **TFRC/TfR1** | BBB shuttles already well established | Mouse purified EC negative in ischemia and GSE225948 day2 pooled EC negative; **receptor abundance not equivalent to productive shuttling** | Use only with real shuttle construct + Fc/affinity design; capillary-depletion, intact cargo bioanalysis, target RNA/protein knockdown, peripheral PK |
| **SLC3A2/CD98hc** | Proven bispecific BBB antibody shuttles in humanized mice/NHPs, distinct uptake and retention compared to TfR (Nat Commun 2023; Cell Chem Biol 2024; follow-up whole-body mapping 2025) | GSE223744 EC roughly stable Slc3a2 RNA; GSE225948 pooled day2 EC p=.743; broad RNA does **not establish** selectivity or CNS accumulation | **Plausible serious comparator** in drug-delivery platform optimization, but not novel receptor identification. Need cargo-access and industrial licensing/partner, stable label, capillary depletion and functional engagement |
| **IGF1R** | Candidate transport handle with peripheral endocrine biology; must establish specific verified shuttle construct and species cross-reactivity | Independent EC down; GSE225948 pooled EC down (p=.0286; q=.133) | Receptor expression or immunofluorescence alone insufficient; require direct validated shuttling platform and safety |
| **LRP1** | Often discussed in targeting and transcytosis; distinct receptor surface pools must be established | Mouse independent EC lower weakly; GSE225948 no detectable overall effect | Need ligand/construct and luminal-to-parenchymal delivery, no blanket assumption of productive shuttle |

## Translation boundary and study design

- A new vascular PET tracer must show measurable binding to an **accessible protein**, normal organ background, dose/timing PK, blocking specificity, radiometabolites, BBB leak/perfusion controls, and tissue-vessel localization. Omics only nominates what to validate.
- A labeled drug or carrier PET study must distinguish parent tracer versus metabolites, carrier versus cargo, capillary adherence versus parenchymal delivery, cell entry versus pharmacodynamic knockdown.
- Existing work makes **CD98hc as a BBB shuttle already established**; a new project must solve a new exposure/target-engagement question, not claim discovery of CD98hc transport.
- If no existing therapeutic/vector is supplied by a scientific/industry partner, do **not** combine complex ASO synthesis + novel shuttle + radiotracer chemistry + rodent stroke surgery into one three-year fellowship.
- For industry career value, prioritize ownership of a quantitative **exposure–distribution–target-engagement** workflow with reproducible analysis over mere multitracer uptake maps.

### Selected verifiable primary literature
- Host AMOT sCD146 PET in myocardial ischemia: https://pmc.ncbi.nlm.nih.gov/articles/PMC8466330/
- CD98hc antibody shuttle in humanized mice and cynomolgus monkeys: https://pubmed.ncbi.nlm.nih.gov/37598178/
- CD98hc bispecific antibody shuttle and tissue-cell targeting: https://pmc.ncbi.nlm.nih.gov/articles/PMC10922565/
- TfR/CD98hc distinct organ and brain cell distribution (2025): https://pmc.ncbi.nlm.nih.gov/articles/PMC11842567/
- BBB-receptor shuttles and validation criteria: https://pmc.ncbi.nlm.nih.gov/articles/PMC12363266/

**Next checks:** protein accessibility/validation; within-original-author-endothelial-subtype GSE225948 sensitivity (Notebook 11 amended and running); systematic tracer/BBB technology literature matrices; industry partner asset confirmation. No tracer or BBB receptor currently earns an unconditional GO.
