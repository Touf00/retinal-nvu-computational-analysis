# IMACS 2026 — proposal screening decision, NOT a dry-lab publication or a human tissue project

Date 13 October 2026. Host candidates: Philippe Garrigue and Benjamin Guillet, CERIMED / Aix-Marseille. This decision memo **supersedes** the earlier suggestion that human pathology procurement, external neuropathologist outreach, PET-ready target validation, or a broad 29-gene study must be completed before drafting an IMACS fellowship. None is a prerequisite to proposal screening. Do not contact tissue banks, laboratories or suppliers for this proposal-screening phase.

## Decision: one project, a staged tracer decision, a genuinely feasible stroke model

**Suggested working title:** *Molecular imaging of neurovascular repair after ischemic stroke: feasibility and prognostic value of existing angiogenesis PET tracers.*

**Core question:** Can an *already developed* angiogenesis PET tracer from the Garrigue/Guillet platform report **delayed neurovascular recovery following focal cerebral ischemia**, independently of injury-induced blood-brain barrier (BBB) leakage and perfusion alterations, and predict vascular/functional recovery?

**Working hypothesis, deliberately testable not asserted:** Early ischemic endothelium and delayed vascular remodeling differ. Although original acute endothelial RNA screens show reductions of AMOT and APLNR in stroke, a **regional/isoform-specific or later postischemic protein signal** might emerge in peri-infarct neovasculature and be imaged. If no specific signal exists in ischemic brain tissue or the tracer cannot distinguish target binding from damaged-barrier accumulation, the hypothesis is rejected quickly and a fallback imaging framework remains feasible.

**Host platform already supported by publications:**
- **[68Ga]Ga-sCD146** targets angiomotin AMOT. Garrigue/Guillet and colleagues used ligand-blocked PET plus ex vivo autoradiography to follow angiogenesis after *mouse hindlimb ischemia* (Theranostics 2018, https://pubmed.ncbi.nlm.nih.gov/30429881/), and subsequent *myocardial infarction* imaging (Cells 2021, https://pubmed.ncbi.nlm.nih.gov/34571954/). **No equivalent rat/mouse ischemic brain PET demonstration was established in reviewed sources.** sCD146 is an AMOT-binding ligand; do not mistake it for an anti-MCAM/CD146 tracer.
- **[68Ga]Ga-AP747** targets APJ/APLNR. Garrigue/Guillet and colleagues radiolabeled and characterized it in Matrigel and mouse *hindlimb ischemia* imaging with pig organ biodistribution (Angiogenesis 2023, https://pubmed.ncbi.nlm.nih.gov/36973482/). No prior cerebral stroke application proved here. APJ receptor is not equivalent to apelin (APLN) ligand expression.
- Garrigue and Guillet also published *transient MCAO in rat* with longitudinal cerebral blood flow, BBB disruption and apoptosis **SPECT/CT** (Int J Stroke 2016, https://pubmed.ncbi.nlm.nih.gov/26763027/). This verifies relevant previous experience, NOT current animal protocol approval.
- CERIMED lists rodent micro-PET/CT, micro-SPECT/CT, PET radiolabeling, autoradiography, optical imaging and histology https://carnotstar.univ-amu.fr/en/technology-platforms/cerimed . **Availability on our project dates, tracer stock, isotope quotas, and small-animal MRI access are NOT confirmed.**
- IMACS provides a 36-month postdoc in interdisciplinary biomedical imaging and requires an international academic secondment and a nonacademic mentor/partner https://euraxess.ec.europa.eu/jobs/458790 . Specific people and placement roles are not yet confirmed.

## Fair comparison of proposal routes

| Route | Fit to existing host technologies | Actual stroke-target evidence | Main cost/feasibility issue | Decision |
|---|---|---|---|---|
| **A. Existing [68Ga]Ga-sCD146 -> AMOT cerebral vascular repair PET** | **High**, host synthesized/tested ligand twice in peripheral and cardiac ischemia | Acute mouse EC Amot RNA decreases in independent datasets; delayed ischemic brain accessible AMOT protein **unproved** | Brain binding, AMOT isoform/luminal availability, ligand plasma residence, BBB leakage | **PRIMARY FEASIBILITY CANDIDATE**; no-go rapidly if autoradiography/brain kinetics fail |
| **B. Existing [68Ga]Ga-AP747 -> APJ cerebral angiogenesis PET** | **High**, host PET ligand, hindlimb ischemia, pig biodistribution | Early mouse EC Aplnr decreases; Apln ligand increases; brain stroke APJ protein dynamics **unproved** | Low human normal EC RNA detection, receptor availability/occupancy, kidney/blood pool, BBB leakage | **PRE-SPECIFIED SUBSTITUTE**, not second full tracer programme |
| **C. ICAM1 or VCAM1 inflammation PET** | Host can image/label, but no host-developed stroke-specific PET radioligand established | Published human stroke vessel protein IHC stronger than AMOT/APJ; not human PET validation | De novo validated binder, chelation, whole-body background, radiotracer pharmacokinetics | **DO NOT MAKE CORE**, use tissue/molecular inflammation comparators if existing reagents |
| **D. Imaging neurovascular flow, BBB change and reperfusion, no new targeting tracer** | Host already published rat tMCAO flow/BBB SPECT; strong CERIMED imaging fit | Well-established stroke endpoints but less molecular target novelty | Need predefined acquisition and reproducible functional endpoint | **GUARANTEED NONMOLECULAR FALLBACK / CO-PRIMARY FUNCTIONAL CONTEXT**, not a claim of angiogenesis imaging |
| E. Broad new GBM or multi-disease target screen | Lab has other disease interests, but no disease-/tracer-matched host and dataset validation for *this* stroke application | Separate project | Dilutes proposal, increases uncertain model/target access | **Not part of this proposal**; revisit only if hosts strongly prefer different disease |

**Important distinction:** Human stroke IHC evidence strength is NOT the same criterion as a practical three-year molecular imaging fellowship built on pre-existing host radiotracers. A more biologically established target can still be a worse project if it requires de novo radiochemistry and major collaborators.

## What the current IMACS archive actually says (screen, not target validation)

The six-animal GSE174574 acute broad EC pseudobulk and independent four-versus-four GSE223744 FACS EC both report:
- **AMOT**: -0.779 and **-1.581** log2CPM mean differences, respectively (acute stroke vs sham, and ischemia vs *naive*). FACS frozen-29-gene BH q ~0.0487. The naive-control confound, frozen-panel testing and small sample sizes limit inference.
- **APLNR (APJ)**: -0.427 and **-1.417**, respectively; FACS q ~0.124. In contrast **APLN (ligand)** rises +2.760 and +4.439 (FACS q ~0.0487). **Ligand upregulation does not establish APJ receptor abundance or PET binding.**
- **ICAM1**: +1.003 and +2.210 (FACS q ~0.0487), supporting its use as an inflammatory endothelial context marker.
- **VCAM1**: +0.246 and +1.118 (FACS q ~0.0487), similarly contextual, not automatically an available tracer.
- Six adult *nonstroke surgical-control* human cerebral FACS EC libraries show median per-sample detected RNA fractions **AMOT ~4.72%, APLNR ~0.12%**. Low capture is a feasibility caution, NOT proof ischemic protein absence.
- GSE225948 pooled stroke EC day14 minus day2 **Amot +0.0185**, nominal p ~0.429 and panel q ~0.600; **Aplnr -0.0117**, p=1.0. They are pooled preparations with batch issues; this does NOT verify delayed angiogenic induction.
- GSE227651 day1,3,7 endothelial profiles suggest only exploratory Amot day7 relative to sham increase while day1/day3 decrease. **One deposited sample per timepoint; no biological replicates, so not evidence of a delayed induction effect.**

**Consequence:** A testable day3-14 protein/PET time course is a *research question*. Do **not** sell AMOT as already upregulated in ischemic brain; do not select APJ solely from APLN RNA. No sequencing analysis must be completed before this grant proposal can be shaped.

Primary local quantitative sources:
- results/target_validation_gate/21_frozen_29_gene_human_mouse_evidence_gates.csv
- results/independent_late/12_GSE225948_day2_EC_exact_permutation.csv
- results/temporal/02_Endothelial_descriptive_ratios.csv
- results/target_validation_gate/22_human_ischemic_vascular_protein_and_tracer_evidence.csv

## A three-aim, 36-month **research project** — imaging first

**Aim 1 (Months 1-6): Fail-fast cerebral tracer feasibility.**
Obtain **host-confirmed** access to one [68Ga]Ga-sCD146 tracer preparation and [68Ga]Ga-AP747 availability **only as conditional backup**; do not presume stock/current protocols. In available *pre-existing experimental stroke or newly approved model* brain samples (sham, acute ischemia and recovery at 3/7/14 days as available), examine AMOT isoforms/EC localization and APJ receptor with ordinary brain-vessel IHC/autoradiography. Test matched radioligand competition on sections and quantify normal brain/vascular background and plasma stability/biodistribution as needed. If no brain-specific tracer signal above binding/permeability controls, **STOP primary tracer and switch once to backup**. A negative feasibility result is a valid milestone, not a three-year failure. Small discovery cohorts do not justify statistical efficacy claims.

**Critical early host decision:** choose the *existing/approved rat transient MCAO* if tracer binding to rat AMOT/APJ is demonstrated and microPET spatial contrast is adequate. If rat ligand cross-reactivity or image resolution invalidates the model, a *mouse* MCAO model may align better with originally validated ligand species, subject to host availability and approvals. Never silently assume both models can run, and do not add unneeded knockout lines.

**Aim 2 (Months 6-20): Longitudinal PET of ischemic neurovascular recovery in ONE rodent model.**
Using the winning tracer, acquire sham-controlled PET at preregistered early/recovery points (e.g. 3,7,14 days), with appropriate dynamic tracer exposure/blood-pool reference, and independent cerebral blood flow and BBB-permeability measurements from host established modalities. Define tissue ROIs from independently mapped lesion/perilesion anatomies; use capillary/vessel immunostaining or autoradiography to distinguish vascular regrowth, mural-cell repair, gliosis and nonspecific injury-leak. Main imaging outcome = specific PET uptake after controlling perfusion, vascular volume and BBB leakage; **not raw lesion SUV**. Correlate signal at a prespecified earlier day with later vascular function, neurological score and relevant tissue-level recovery readout (e.g. day28). Analysis unit = **animal**, repeated scans nested within animal. Estimate powered n after pilot tracer variance; do not invent an n now.

**Aim 3 (Months 20-30): Does molecular PET add functional prediction or treatment-monitoring value?**
Use an *already established within-host* animal intervention model: transient reperfusion versus sustained occlusion if ethically feasible, or **one** previously validated host-available therapy (pharmacological, oligonucleotide or vascular intervention) with a fixed mechanism and vehicle/sham controls. Do not promise access to a therapy not named/verified by hosts. Compare whether early PET predicts later perfusion, vascular remodeling and behaviour **after accounting for initial lesion size and flow**. For a drug study, track receptor occupancy/competition and physiological effects of therapeutic agent so ligand displacement is not mistaken for target downregulation.

**Months 30-36:** Analyze repeated-measure scans with animal-level mixed effects, quantitative kinetic/partial-volume and permeability correction, robust image registration, annotated reproducible workflows, translational decision report, articles/abstracts if justified. Integrate fixed *3-month international academic secondment* (quantitative PET kinetics or imaging AI rather than repeating biology) and a **named nonacademic mentor** (radiopharmaceutical industrial production, imaging CRO or clinical imaging setting). These roles require confirmation in the actual IMACS application.

## Essential feasibility gate and reviewer objections

**DANGER 1: AMOT/APJ early down, not up.**
Answer: peripheral ischemia literature provides tracer platform, **not a CNS hypothesis proven in advance**. Late-stage tissue protein and brain PET are the actual test. Abandon if absent; do not force a new target.

**DANGER 2: Brain tracer is simply BBB leakage.**
Answer: matched dynamic kinetics, blood-input/blood-pool correction, BBB injury measurement, radioligand-specific blocking/autoradiography and comparison with nonbinding probe when possible. If leakage cannot be separated from specific binding, NO claimed molecular stroke tracer.

**DANGER 3: MR/PET small brain and species limits.**
Answer: micro-PET hardware is available at CERIMED, but current scheduling and resolution must be checked. Choose rat versus mouse only after binding + physics pilot; co-register with reliable brain anatomy, and use ex vivo validation. Do not presume rodent MRI from a listed 3T large-animal scanner.

**DANGER 4: Not enough time, too many models/targets.**
Answer: maximum ONE disease, ONE animal model, ONE primary tracer; backup chosen only at six-month go/no-go; ONE functional/therapeutic perturbation. Three actual work packages, not nine candidate-specific clinical studies or new human tissue collection.

**DANGER 5: Therapeutic meaning.**
Answer: longitudinal molecular imaging seeks a *predictive* readout of vascular recovery/reperfusion and treatment response, not a new medicine. A therapy intervention is conditional on a concrete host-supported choice and is not falsely promised in application.

**DANGER 6: Clinical translation.**
Answer: proof of concept is preclinical. Existing human normal EC and postmortem bulk stroke RNA contextualize relevance but cannot validate target accessibility. No human autopsy/brain bank outreach is needed to submit or screen this proposal.

## What is needed from Garrigue/Guillet before finalizing (FOUR factual checks, not new collaborations)

1. Is **[68Ga]Ga-sCD146** reproducibly synthesizable at CERIMED now, or is **[68Ga]Ga-AP747** the better currently available ligand? What are recent yield/stability and animal PET availability? Paper does not guarantee stock.
2. Which **single** ischemic stroke model is routinely operational (rat transient MCAO vs mouse MCAO), and are sections, tissue histology and perfusion/BBB readouts available for a quick brain-specific binding/contrast pilot?
3. What is the **one concrete** reperfusion/therapeutic perturbation and recovery endpoint their current collaborators can implement, without setting up a new disease model?
4. Which **international academic** and **non-academic** IMACS partner are truly available for imaging quantification/radiopharmacy/industry skills? IMACS explicitly requires both; never invent placements.

These can be resolved in a short meeting or email with the already engaged hosts. **No external neuropathologist, biobank permission, human tissue transfer, company outreach or new lab procurement is a proposal-screening requirement.**

## User career value

Hands-on longitudinal **PET/SPECT tracer quantification**, tracer chemistry and ligand specificity evaluation, **stroke wet-lab neurovascular physiology**, image registration/kinetic modeling, preclinical theranostic decision gates, animal-level study design, translation-to-industry imaging biomarker workflow, and collaboration with quantitative imaging scientists. This deliberately delivers portable pharma/biotech skills while preserving the user's experimental neuroscience identity; do not market this as primarily a dry-lab paper.

## Bottom line

**Proceed with a one-tracer neurovascular repair imaging proposal, anchored in the host's existing [68Ga]Ga-sCD146 and conditional [68Ga]Ga-AP747 PET tracers, and verified rat cerebral ischemia imaging experience.** The six-month brain feasibility gate is not a preliminary human-stroke protein study. Reject brain targeting early if inaccessible, and maintain a meaningful imaging/perfusion outcome. This is a strong *application hypothesis*, not a guaranteed successful tracer result.

This screen is the CURRENT proposal direction; prior files about human tissue procurement are **optional hypothetical later validation ideas**, NOT tasks to perform now. All edits on isolated branch imacs-stroke-pet-screen, no merge to main.
