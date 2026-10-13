# START HERE — IMACS fellowship application control document

**Last consolidated:** 14 October 2026. **Use this as the single index; older conversations and exploratory candidate lists are not separate commitments.**

## 1. Deadline and official programme

- IMACS Marseille/Aix-Marseille postdoctoral application deadline: **30 November 2026, 23:59 Central European local time** (official EURAXESS posting displays `Europe/Vienna`, same civil clock as Paris).
- Fellowship: 36 months; advertised appointment start 1 September 2027.
- Host: Philippe Garrigue and Benjamin Guillet at Aix-Marseille/CERIMED, **subject to host commitments and final nomination**.
- Required: foreign academic secondment **3 months**, and one nonacademic partner/mentor.
- Sources: https://euraxess.ec.europa.eu/jobs/458790 ; https://institut-imaging.univ-amu.fr/en/imacs-cofund-programme

**Verify exact attachment templates, page inclusions, letter signatories, formatting, portal submission details and allowed budgets against the current official IMACS applicant guide.** Official host portal was not retrievable for a line-by-line check on 14 October.

## 2. Current document checklist

From the **previous IMACS call-guide review**, pending fresh template validation:

- Research and Training Proposal — **up to 10 pages**.
- CV — **up to 5 pages**.
- Data Management Plan — **1–2 pages**.
- PhD certificate / proof of research eligibility; ID/passport; ethics self-assessment if applicable.
- Commitment letters: supervisor; host laboratory; international academic secondment institution; nonacademic mentor/partner.
- Online submission; not simply an informal email to the host.

The EURAXESS call independently verifies foreign partner, nonacademic partner, three-month academic secondment, 36-month position, and deadline. **Do not use the ordinary European Commission MSCA COFUND organisation-level grant application template** in place of this fellowship's institution-specific applicant template.

## 3. Scientific direction to discuss with Garrigue / Guillet NOW

**Provisional title:** *Molecular PET imaging of neurovascular repair after ischemic stroke: feasibility and predictive value of existing angiogenesis tracers.*

**Question:** Can an already available host tracer distinguish delayed vascular repair after cerebral ischemia from disrupted flow, blood-brain barrier leakage, and nonspecific tracer accumulation, and improve prediction of subsequent recovery?

**Primary feasibility candidate:** [68Ga]Ga-sCD146 ligand targeting **angiomotin (AMOT)**, not MCAM/CD146. **Conditional alternative:** [68Ga]Ga-AP747 targeting **APJ/APLNR**. This ordering is provisional. Host chemistry availability, validated brain luminal accessibility, and species cross-reactivity decide.

**Main warning from computation:** Independent mouse acute endothelial cohorts show **AMOT and APLNR RNA reductions after MCAO**. Neither is established as a PET-accessible cerebral post-stroke protein. Do not promise increased APJ or AMOT PET signal.

**Core experimental structure:**
- WP1 / Months 1–6 — species/brain tracer specificity, protein localization, ex-vivo blocking/autoradiography, BBB confound controls. Predeclared GO/NO-GO. Choose ONE tracer or fall back to physiologic perfusion/BBB imaging.
- WP2 / Months 6–20 — longitudinal PET with perfusion and BBB measurements in ONE active host stroke model, lesion-defined ROIs and blinded tissue validation; animal is unit of inference.
- WP3 / Months 20–30 — test whether early imaging predicts later vascular/neurologic recovery, beyond baseline lesion injury/CBF; one host-supported treatment/reperfusion perturbation only if feasible.
- Months 30–36 — kinetic/statistical validation, reports, writing, handoff, reproducible image-analysis pipelines. Include three-month foreign academic secondment and real industrial mentorship within the timetable.

**Do not conflate a useful proposal with a completed dry-lab publication.** No human tissue bank outreach or extra clinical PI permission needed to choose the current fellowship hypothesis.

More detailed concept: [one-page host discussion](IMACS_2026_ONE_PAGE_PET_STROKE_CONCEPT.md).
Decision rationale: [red-team proposal screen](IMACS_2026_PROPOSAL_SCREEN_AND_PROJECT_DECISION.md).

## 4. Required partner decisions — contact hosts promptly

Ask Garrigue/Guillet to confirm:
1. What tracer can be reliably produced **now**: sCD146/AMOT and/or AP747/APJ? Any unpublished brain-stroke data?
2. Which **one** stroke model is currently active (rat tMCAO, mouse, or other), which perfusion/BBB readouts are schedulable, and what relevant tissue is already available?
3. One feasible disease/recovery/therapy comparison. Avoid inventing an inaccessible therapeutic solely to make the project sound industry-oriented.
4. Which international academic collaborator can provide three months of **quantitative PET kinetics, translational imaging or multimodal image analysis**; who is the supervisor?
5. Which nonacademic organisation can provide named **radiopharmaceutical/biomarker or drug-development** mentorship? AAA/Novartis is one option; not the only permissible partner. Choose confirmed scientific contribution, not a logo.
6. Who signs host and supervisor commitment letters, and which IMACS templates are required.

**Do not wait until all datasets are fully resolved before this host contact.**

## 5. Recommended 10-page proposal architecture — NOT a quoted official template

1. Title/abstract, unmet need, central hypothesis and objectives — 0.7 p.
2. Prior art and specific imaging gap, host capability, preliminary research — 1.3 p.
3. Three work packages / controls / animal-level design / quantitative imaging — 4.0 p.
4. Risk, alternatives, validation, technical go/no-go criteria — 0.8 p.
5. Applicant training and supervision / secondment / industry mentorship — 1.1 p.
6. Novelty, impact, translation, career development, dissemination — 0.8 p.
7. Milestones, deliverables and 36-month **Gantt chart** — 0.7 p.
8. References — 0.6 p.

Sum: 10 pages. Reallocate once the official template specifies figure/reference inclusion and headings. A Gantt chart is **strongly recommended for reviewer clarity**, but its status as a formal mandatory attachment is not independently established.

### Simple 36-month timeline
| Work package | M1–6 | M7–12 | M13–18 | M19–24 | M25–30 | M31–36 |
|---|---|---|---|---|---|---|
| Tracer feasibility + no-go gate | X | | | | | |
| Longitudinal PET/perfusion/BBB | | X | X | X | | |
| Recovery/intervention prediction | | | | X | X | |
| Academic secondment | | X or later | X or later | | | |
| Nonacademic mentor engagement | X | X | X | X | X | X |
| Analysis, reproducibility, dissemination | X | X | X | X | X | X |

## 6. Canonical computation and limitations

Code and successful processed results are committed on this **isolated branch inside an EXISTING public retinal repository**, not a dedicated IMACS repo.

- [Recovery instructions](RECOVERY_AND_STATUS.md) — clone commands, eight original notebooks plus further extensions, historical runs and permanent evidence tables.
- [Critical scientific audit](results/SCIENTIFIC_AUDIT.md) — flags false replication and unvalidated annotation.
- [Results directory](results/) — independent human/mouse, rat, spatial, atlas and integration.
- [Human vascular protein / tracer audit](results/target_validation_gate/22_human_ischemic_vascular_protein_and_tracer_evidence.csv).

**Crucial limitations:**
- n3+3 primary MCAO: exact two-sided balanced label-permutation minimum p=.10.
- Unreplicated day1/3/7 mouse trajectory; no animal-level inference.
- Rat hemisphere preparations are pooled, and some endothelial labels are very sparse.
- Visium mixed-cell spots lack validated author infarct region annotations; do not claim core-versus-penumbra specificity.
- SCP3078 complete 5GB atlas was inaccessible in the automated environment; only metadata audited.
- Mouse FACS EC independent n4 vs n4 uses **naive rather than sham** controls.
- Human bulk stroke autopsy signatures do NOT prove which molecule is present on the **luminal surface** of human ischemic vessels.
- RNA signals do NOT prove protein expression, successful PET tracer binding, BBB transfer or target engagement.

**No need to repeat analyses simply because the prior conversation ended.** A 'green' workflow shows code execution, not scientific validity.

## 7. Next actions and deadlines

- **Within the next few days:** send the [one-page host discussion concept](IMACS_2026_ONE_PAGE_PET_STROKE_CONCEPT.md) and request 20–30-minute decision call. Ask only the host asset questions.
- **Before end of October:** settle primary/backup tracer logic, one stroke model, partner choices, confirmed signatories and official template checklist. Draft the 10-page research + training plan.
- **First half November:** circulate full proposal for host/partner comments; create a Gantt figure, controls/risk figure, DMP, CV and partner letters.
- **Before 30 November:** quality-control every annex, upload early, verify receipt. Do not target last-hour submission.

**Clinical/translational value to the user:** molecular imaging + experimental neurovascular biology + kinetic modeling, quantitative biomarker/PK-PD decision-making and real nonacademic mentoring. The fellowship must deliver both industry-portable methods and a personally owned scientific question.

## References / source priority

1. Official IMACS Marseille 2026 EURAXESS call: https://euraxess.ec.europa.eu/jobs/458790
2. IMACS institutional programme and templates: https://institut-imaging.univ-amu.fr/en/imacs-cofund-programme
3. Specific internal audited results and primary papers linked inside proposal-screen documents.

**If a detail in this index conflicts with the officially downloaded IMACS guide/template, the official guide controls.**