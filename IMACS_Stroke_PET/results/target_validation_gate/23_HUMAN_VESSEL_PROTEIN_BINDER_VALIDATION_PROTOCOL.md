# IMACS human ischemic vessel protein and binder validation pilot

Status: REPRODUCIBLE PLAN AND ANALYSIS CODE PREPARED; NO PATIENT ASSAY PERFORMED. Prepared 13 October 2026. Targets ICAM1/CD54 and VCAM1/CD106. Other target classes stay in the existing 9-target hypothesis matrix, not simultaneous radiotracer development.

## Scientific premise

Human infarct tissue immunohistochemistry (IHC) documents ICAM1-positive microvessels in 11 patients (PMID 8790029) and infarct VCAM1 in endothelial cells AND astrocytes in six patients (PMID 10456619). But IHC, particularly formalin-fixed sections with antigen retrieval, cannot demonstrate whether a bloodborne ligand binds the extracellular luminal epitope of a living ischemic brain vessel.

Gate A = endothelial-localized HUMAN ischemic brain protein, compared within donor and by lesion region. Gate B = independently demonstrated specific binding to LIVE, nonpermeabilized, blood-facing human BBB endothelium under controlled flow. PET progression requires BOTH; protein alone, RNA, and animal MRI do not pass PET validation.

## Gate A: tissue-access and blinded human pathology

### Access and approvals

- FIRST REQUEST: BRAIN UK https://brain-uk.org/home/informal-enquiry/ ; brainuk@soton.ac.uk. It explicitly includes cerebrovascular disease among tissue types, is open to overseas institutional applicants, and primarily supplies archived FFPE. Its database is not preconfirmed to hold the needed matched infarct region blocks. A bona fide institutional PI, ethical and local research approvals, sample-transfer agreement and cost coverage are needed. Samples may take months even after approval.
- LOCAL MARSEILLE: CRB AP-HM https://fr.ap-hm.fr/site/crb-centre-de-ressources-biologiques ; CRB.Timone.Marseille@ap-hm.fr. ASK whether a neuropathology collaboration can locate autopsy infarct cerebral FFPE blocks. Its published sample collections do NOT establish that it has acute ischemic cerebral samples.
- Backup: Netherlands Brain Bank via e-NBB; stroke cohort not confirmed. Neuro-CEB primarily neurodegeneration and controls, so not an acute-stroke sourcing guarantee.
- No transfer, use, or public disclosure of donor clinical information without appropriate permissions. Store ONLY de-identified vessel counts and summary tables in GitHub. Raw patient images and linked clinical identifiers stay in approved institutional systems.

### Minimum target sample request

For pilot feasibility seek about 8-12 independent adult acute ischemic stroke autopsy patients (documented onset-to-death 12-72 h, including uncertainty) with neuropathologist-mapped infarct core, border/peri-infarct and anatomically matched contralateral brain blocks wherever actually available. Seek 5-8 nonstroke brain control cases matched for age, brain region and fixation as a separate exploratory comparison; no claim this n is statistically powered. Control and stroke tissue must share comparable brain region/vascular segments and imaging batch. Subacute cases over 72 hours and hemorrhagic transformation are labeled and analyzed separately.

Record anonymized donor ID, diagnosis, age/sex, clinical interval, reperfusion/thrombectomy, sepsis/systemic inflammation, hemorrhagic conversion, stroke territory, onset uncertainty, postmortem interval, brain region, block, hemisphere, fixation duration, sample preparation and staining batch. Pathology team, not the ICAM1/VCAM1 fluorescence, must establish core/border/contralateral identity. Never infer unavailable contralateral tissue, a precise core from anatomy alone, or cross-donor pairs.

### IHC reagents and compartment identification

Qualified candidate human FFPE clones, NOT human luminal binders:

- ICAM1: EPR4776, Abcam ab109361, human IHC-P and manufacturer KO controls. Verify human ischemic brain assay on a known positive tissue and independent antibody/antigen loss control; vendor example is not brain validation.
- VCAM1: EPR5047, Abcam ab134047, human IHC-P and manufacturer KO controls. Test positive human tonsil/cytokine activated endothelial pellet and clean negative controls.

Separate primary-target stains/panels with endothelial ERG nuclear and CD31/PECAM1 or CDH5 membrane markers, collagen IV/laminin vascular basement membrane. A second panel uses PDGFRB/NG2 pericytes, GFAP astrocytes and CD45 leukocytes. Optional albumin/IgG/fibrinogen for BBB damage, glycophorin A for intravascular erythrocytes. A target-positive leukocyte stuck to endothelium must NOT be counted as endothelial protein. Tissue autofluorescence/hemosiderin, necrotic core and loss of endothelial area must be documented. Spectral bleedthrough and species/antibody channel incompatibility must be optimized, potentially on registered serial sections.

Predefine systematic image sampling by region (no positive-field selection), equal acquisition settings within batches, randomized staining order, no-primary/secondary-only/isotype and positive-control specimens; calibrate with three or more independent negative control fields for each target and batch. Record the target-protein cutoff BEFORE clinical region unblinding and never optimize cutoff to maximize stroke contrast.

Experienced blinded assessors independently classify vessels (capillary, postcapillary venule, arteriole or indeterminate) and verify mask quality and true pathology region. Vessel ROI = actual endothelial membrane/ring signal, not entire perivascular space. Indeterminate or small masks are tracked and excluded according to protocol. Consider inter-rater reannotation of 15% random ROIs as qualitative reliability check.

### Primary and secondary endpoints

Primary: in ischemic patients whose time-to-death is <=72 hours, **donor-level difference in proportion of ICAM1- or VCAM1-positive endothelial capillaries**, peri-infarct minus matched contralateral brain. Two targets are co-primary exploratory and the only two comparisons in BH if both have >=6 paired donors. A patient must have >=10 adequately segmented capillaries in EACH paired region for a target; otherwise the patient is excluded from that target's paired test, with missingness reported. Capillary endpoint addresses BBB imaging; postcapillary venules must be evaluated as a distinct prespecified SECONDARY segment because adhesion markers may favor those vessels.

Other descriptive endpoints: infarct core separately; endothelial vessel density per region and mean background-corrected endothelial target signal; fraction of target signal in nonendothelial GFAP+/CD45+/mural compartments; contrast by lesion time and hemorrhagic conversion. Never claim vessel-level repeated images are independent patients. No statistical human protein result exists yet.

Analysis script: IMACS_Stroke_PET/scripts/analyze_human_ischemic_vessel_IHC.py. It requires prelocked staining thresholds and a deidentified vessel-level CSV. It uses donor-level exact paired signflips (n>=6, otherwise descriptive only) and donor-resampling bootstrap intervals. For n=6 the theoretical minimum attainable two-sided exact p=0.03125; this pilot has no guarantee of FDR discoveries or statistical power.

CSV fields for vessels:
donor_id,block_id,slide_id,vessel_id,target,assay,stain_batch,clinical_group,poststroke_hours,region,segment,endothelial_mask_pixels,target_mean_intensity,local_background_intensity,mask_qc_pass,pathology_region_verified,acquisition_blinded

Locked calibration CSV:
target,stain_batch,corrected_intensity_cutoff,n_independent_negative_control_fields,positive_control_pass,threshold_locked_before_unblinding

Controlled values: assay=FFPE_IHC; clinical_group AIS or nonstroke_control; regions peri_infarct, infarct_core, contralateral, nonstroke_control; segments capillary, postcapillary_venule, arteriole, indeterminate. Boolean fields exactly True/False strings. One vessel ROI per antigen per unique slide key; all source IDs deidentified.

Run using real external data only after tissue access:
    python IMACS_Stroke_PET/scripts/analyze_human_ischemic_vessel_IHC.py --vessels /approved/path/vessels.csv --thresholds /approved/path/locked_control.csv --out /approved/path/results/

Synthetic self-tests only, NOT clinical results, are implemented in GitHub Actions workflow imacs_human_vessel_ihc_safety.yml.

## Gate B: separate human BBB luminal antigen binding

Gate B cannot be run on FFPE. Use human brain microvascular endothelial cells (primary HBMECs or credible iPSC-derived BBB endothelium) grown in a lumenized perfused 3D BBB chip with collagen/laminin support and, when possible, pericyte/astrocyte co-culture. Verify endothelial identity (PECAM1, CLDN5, CDH5), viability, continuous lumen and functional tracer permeability BEFORE target labeling. Peer-reviewed model: https://www.nature.com/articles/s41596-021-00635-w ; human BBB TNF induced ICAM1/VCAM1 under flow: https://pmc.ncbi.nlm.nih.gov/articles/PMC9610437/

Use at least three independent biological cell donors/differentiation batches for a feasibility study, with more technical chips per independent donor; chips are not independent humans. Compare a well-defined control, defined cytokine activation and, if mechanistically justified, an additional OGD/reoxygenation condition with independently monitored viability and barrier permeability. TNF activation models inflammation, not clinical stroke tissue.

Live-surface candidate reagents need a DIFFERENT validation from FFPE clones: anti-human ICAM1 HA58 (BD 555510) and anti-human VCAM1 51-10C9 (BD 561679), supported for human flow/surface binding. Commercial labeled formulations may contain azide/BSA and are NOT directly suitable for living perfused cell experiments. Acquire a confirmed azide-free, low-endotoxin, non-toxic format, confirm direct binding to the accessible extracellular epitope, lot, blocking effects and fluorescent probe compatibility. A known antibody clone used in flow does NOT automatically qualify as a PET radiotracer.

Perfuse defined labeled targeting ligand under known shear, matched dose and time, and measure retained luminal probe after washing. Compare matched-size and fluorophore/isotope nonbinding probe, excess unlabeled antigen-specific binder competition, immunogen/antigen knockdown where feasible, and endothelial viability. Quantify binding with plasma/perfusate input exposure, and parallel 10/70-kDa dextran (or size-matched tracer) permeability. Localization must distinguish apical lumen from basolateral receptor or leaked tracer in surrounding matrix. Antibody inhibition can affect receptor internalization or flow; those are measured, not assumed absent.

Outcome GO requires reproducible antigen-dependent live luminal binding over a matched nonbinder *under barrier-competent conditions*, distinguishable from simple blood-pool/permeability enrichment, and expression observed on ischemic human vessels in Gate A. If staining occurs only after fixation/permeabilization, signal tracks permeability, a nonbinding probe behaves similarly, or antigen-loss/competition does not reduce signal: no specific accessible cerebral tracer claim.

## Gate C: PET only if A+B support development

Ask radiopharmacy to review existing binders/fragment affinities and a PET isotope whose half-life matches binder kinetics (fast small-fragment vs long-circulating antibodies). Measure radiochemical purity, molar activity, parent radio-metabolites, target competition, organ sinks, in-vivo brain vessel/blood distribution. For stroke model proof randomize at the animal level, compare sham and ischemia with dynamic PET, blood-pool and perfusion references, paired leakage measurement, matched nonbinding ligand and postmortem microvessel autoradiography. Ischemic brain SUV alone never proves accessible antigen.

A BBB therapeutic shuttle is a different endpoint: normal human brain receptor protein (e.g. CD98hc or TFRC) is not stroke selectivity; actual delivery requires intact cargo beyond capillaries, capillary depletion and cell-specific pharmacodynamics.

## Feasibility decisions and blockers

- If institution/biobank cannot provide pathologist-labelled acute stroke with matched regions, only descriptive protein or separate positive-control studies are possible; DO NOT claim donor-paired human stroke findings.
- If endothelial staining reproducible but no live luminal access, publish only localization; STOP PET vascular-target claim.
- If live luminal target binding positive but BBB disruption and matched nonbinding probe explain signal, STOP until confounds resolved.
- If patient count below six paired complete cases, report patient-level descriptive data and donor bootstrap sensitivity, not underpowered p-values as biological proof.
- Tissue availability and ethics are currently NOT confirmed; no reagents ordered, no biohazardous material transferred, no human assay executed.

## Public source anchors

Human ischemic protein: https://pubmed.ncbi.nlm.nih.gov/8790029/ ; https://pubmed.ncbi.nlm.nih.gov/10456619/
Archived tissue request: https://brain-uk.org/home/informal-enquiry/ ; https://brain-uk.org/home/faqs/
APHM Marseille: https://fr.ap-hm.fr/site/crb-centre-de-ressources-biologiques
IHC reagents: https://www.abcam.com/en-us/products/primary-antibodies/icam1-antibody-epr4776-ab109361 ; https://www.abcam.com/en-us/products/primary-antibodies/vcam1-antibody-epr5047-ab134047
Separate live-surface reagents: https://www.bdbiosciences.com/en-se/products/reagents/flow-cytometry-reagents/research-reagents/single-color-antibodies-ruo/purified-mouse-anti-human-cd54.555510 ; https://www.bdbiosciences.com/en-us/products/reagents/flow-cytometry-reagents/research-reagents/single-color-antibodies-ruo/pe-mouse-anti-human-cd106.561679

Source audit files: 23_human_stroke_tissue_access_routes.csv and 23_antibody_qualification_shortlist.csv under results/target_validation_gate/. Keep all work on imacs-stroke-pet-screen; do not merge to retinal main.
