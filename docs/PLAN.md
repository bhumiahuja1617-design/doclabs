# DockLab Implementation Plan & Roadmap

**Author**: Senior Full-Stack Engineer & Computational Chemist  
**Project**: DockLab (AutoDock Vina Docking, Automated Structure QC, and Multi-Seed Reliability Engine)  
**Specification**: DockLab Master Prompt v2 with Biotech Depth Layer (Section 19)

---

## 1. Executive Summary & Core Identity

DockLab is a high-rigor computational biochemistry web application and science engine designed to execute protein-ligand docking with AutoDock Vina, while addressing the fundamental limitation of standard docking tools: **telling the researcher whether the computational result can actually be trusted.**

DockLab is defined by two signature features:
1. **Pre-docking Structure & Parameter QC Engine**: Rigorously inspects receptor coordinates (missing heavy atoms, gaps, bridging waters, cofactors/metals), ligand stereochemistry, flexibility, charge states, and search box boundaries before a single docking step is computed, emitting actionable PASS / WARN / FAIL diagnostics.
2. **Multi-Seed Docking Reliability Engine**: Re-docks top-ranked candidate ligands across multiple independent pseudorandom seeds ($N=5 \dots 10$), evaluating score variance ($\sigma \le 0.3 \text{ kcal/mol}$), pose clustering RMSD consistency ($\ge 80\%$), interaction conservation, and Spearman rank stability to distinguish genuinely robust poses from stochastic search artifacts.

All computational outputs are explicitly framed as **computational repeatability and physical interaction estimations**, never fabricated as "binding probabilities" or "guaranteed real-world efficacy."

---

## 2. Phase-by-Phase Task Breakdown & Roadmap

### Phase 0: Foundations & Environment Verification
- **Goal**: Establish the repository infrastructure, container/environment definitions, automated preflight validation, planning, and foundational biological documentation.
- **Tasks**:
  1. Initialize repository files: `.gitignore`, `.env.example`, pre-commit hooks, pinned `environment.yml`, `backend/Dockerfile`, `frontend/Dockerfile`, `docker-compose.yml` (API, worker, redis, postgres, frontend).
  2. Implement `scripts/preflight.py` verifying critical scientific dependencies: AutoDock Vina, RDKit, Open Babel, Meeko, fpocket, PDBFixer, OpenMM, Biopython, ProLIF, spyrmsd, WeasyPrint, FastAPI, Celery, SQLAlchemy.
  3. Formulate `docs/DECISIONS.md`, `docs/BUILD-LOG.md`, `docs/GLOSSARY.md`.
  4. Write initial learning guide chapters:
     - `docs/learn/00-big-picture.md`
     - `docs/learn/01-environment-and-tools.md`
     - `docs/learn/13-protein-structure-and-binding-basics.md`
     - `docs/learn/14-hiv-protease-and-drug-resistance.md`
- **Acceptance Criteria**: Preflight script runs clean, imports all libraries, executes binary CLI checks (`vina --help`, `fpocket -h`, `obabel -V`), prints verified versions, and exits 0.
- **Approval Gate**: **YES** (stop and present `PLAN.md` and preflight output to project owner).

---

### Phase 1: Science Core CLI & Biological Annotations
- **Goal**: Implement zero-web, standalone scientific calculation modules and verify with CLI execution.
- **Tasks**:
  1. `app/services/protein_service.py`: RCSB PDB fetcher/cache, Biopython parser, chain filtering, alternate location cleanup, HETATM extraction, PDBFixer heavy atom completion and protonation at pH 7.4, Meeko receptor prep to PDBQT (with Open Babel fallback).
  2. `app/services/ligand_service.py`: SMILES parsing, largest fragment isolation, neutralization, pH 7.4 protonation, ETKDGv3 3D embedding with fixed seed, MMFF94 energy minimization, Meeko PDBQT export, drug-likeness computation (Lipinski Rule of 5, Veber rules, PAINS/Brenk alert filtering, ligand efficiency).
  3. `app/services/box_service.py`: Bounding box generation centered on crystal ligand centroid + 10 Å padding (clamped 15–30 Å), residue-based box, fpocket top pocket integration, and user-specified manual bounds.
  4. `app/services/docking_service.py`: AutoDock Vina Python wrapper supporting exhaustiveness presets (Fast=4, Balanced=8, Thorough=32), multi-pose generation, PDBQT and SDF pose export.
  5. `app/services/validation_service.py`: Redocking pipeline using crystal ligand coordinates, symmetry-aware heavy-atom RMSD calculation, and ProLIF key-residue contact recovery.
  6. Target cards implementation: `examples/1hsg/target.json` with verified UniProt, function, catalytic triad/dyad annotations, and known inhibitors.
  7. CLI entrypoint: `scripts/cli_dock.py` for end-to-end headless validation.
  8. Unit & regression tests: `backend/tests/test_science_core.py`.
  9. Learning chapters: `02`, `03`, `04`, `05`, `06`, `15`, `16`, `17`, `18`.
- **Acceptance Criteria**: `cli_dock.py` redocks 1HSG crystal ligand (indinavir/XK263) achieving heavy-atom RMSD $< 2.0\text{ \AA}$ and key-residue recovery $\ge 70\%$.
- **Approval Gate**: **YES** (show real redocking numbers, RMSD, and key-residue recovery).

---

### Phase 2: Pre-Docking Quality Control (QC) Engine
- **Goal**: Implement automated diagnostic inspections across protein, ligands, and bounding box prior to docking.
- **Tasks**:
  1. `app/services/qc_service.py`:
     - Missing residues check (REMARK 465 vs ATOM sequence; gap distance to box).
     - Missing heavy atoms check (REMARK 470, PDBFixer missing list within 8 Å of box).
     - Chain & assembly check (REMARK 350 vs selected chains, contact loss).
     - Bridging waters check (waters within 3.5 Å of ligand and $\ge 2$ protein atoms).
     - Cofactors & metal ions check (HETATM within 6 Å of pocket).
     - Ligand stereochemistry check (`FindMolChiralCenters` for undefined chiral centers).
     - Ligand flexibility and charge sanity (MW, rotatable bonds, formal charge).
     - Box sanity (dimensions, empty volume fraction, pocket overlap, crystal inclusion).
     - Pocket protonation notes (ionizable residues HIS/ASP/GLU/LYS).
  2. Implement typed schemas (`QCReport`, `QCCheck`, severity `pass|warn|fail`).
  3. Fixture testing: run against clean 1HSG (clean baseline) and flawed test fixtures (truncated pocket residue, unassigned racemic ligand).
  4. Integrate `--qc` flag into `cli_dock.py`.
  5. Learning chapter: `07-qc-engine.md`.
- **Acceptance Criteria**: Flawed inputs trigger expected warnings and fails; clean 1HSG triggers no false failures; acknowledgement logic enforced.
- **Approval Gate**: No.

---

### Phase 3: Docking Reliability Check Engine
- **Goal**: Implement multi-seed docking repeatability analysis for top-K candidates.
- **Tasks**:
  1. `app/services/reliability_service.py`:
     - Re-dock top-K ligands ($K=10$) across $N$ seeds ($N=5 \dots 10$).
     - Compute score mean and sample standard deviation ($\sigma$).
     - Calculate pairwise symmetry-aware heavy-atom RMSD across top poses, cluster at $2.0\text{ \AA}$ threshold, and determine largest cluster fraction.
     - Compute rank variation across seeds and mean Spearman rank correlation coefficient ($\rho$).
     - Calculate ProLIF interaction fingerprint conservation across consensus poses.
     - Assign category labels (`High`, `Moderate`, `Low`) with deterministic plain-English explanations.
  2. Unit tests with synthetic pose coordinates and ranking permutations.
  3. CLI flag `--reliability` in `cli_dock.py`.
  4. Learning chapter: `08-reliability-check.md`.
- **Acceptance Criteria**: Executes on 1HSG example; produces valid statistics, cluster metrics, and clear diagnostic explanations.
- **Approval Gate**: No.

---

### Phase 4: Analysis & Two-Target Benchmark
- **Goal**: Complete interaction analytics, virtual screening metrics, and execute two-target benchmark.
- **Tasks**:
  1. `app/services/analysis_service.py`:
     - Contact detection via ProLIF (H-bonds, hydrophobic, $\pi$-stacking, salt bridges).
     - Virtual screening evaluation: ROC curve, AUC, Enrichment Factor (EF1%, EF5%, EF10%).
     - Property bias warnings (MW distribution comparison between actives and decoys).
  2. Prepare second target: EGFR kinase domain (e.g., PDB 1M17 / 4HJO / erlotinib complex) with curated actives and property-matched decoys.
  3. Create `benchmark/run_benchmark.py`: runs redocking, virtual screening, and reliability sanity check on both 1HSG and EGFR.
  4. Generate `benchmark/RESULTS.md` with measured numbers.
  5. Learning chapters: `09-roc-auc-enrichment.md`, `12-benchmark.md`.
- **Acceptance Criteria**: Automated benchmark script completes end-to-end on both targets without hardcoded outputs; generates ROC-AUC and EF5% with real data.
- **Approval Gate**: No.

---

### Phase 4.5: Biological Case Study (HIV-1 Protease Drug Resistance)
- **Goal**: Execute the primary scientific research deliverable defined in Section 19.1.
- **Tasks**:
  1. Curate verified approved HIV-1 protease inhibitors (darunavir, saquinavir, indinavir, ritonavir, nelfinavir, lopinavir, atazanavir) with PubChem CIDs.
  2. Source wild-type (1HSG) and verified experimentally solved crystal mutant structures (e.g., V82A, I50V, I84V, L90M, D30N) from RCSB PDB with resolution $\le 2.5\text{ \AA}$, logging to `docs/case_study/structures.csv`.
  3. Curate published fold-change resistance data from Stanford HIV Drug Resistance Database into `docs/case_study/reference_data.csv` with source URLs.
  4. Implement `benchmark/run_case_study.py`:
     - Docks all inhibitors against wild-type and mutant variants.
     - Calculates $\Delta\text{score}$, contact changes at mutated residues, H-bond count deltas, steric clashes, key-residue recovery, and reliability labels.
     - Computes Spearman rank correlation between $\Delta\text{score}$ and $\log(\text{fold-change})$.
  5. Write `docs/CASE_STUDY.md` (structured paper format) and learning chapter `19-case-study-walkthrough.md`.
- **Acceptance Criteria**: Script runs fully, outputs tables and figures, notes all citations and URLs, reports real correlations honestly (including negative/weak findings).
- **Approval Gate**: **YES** (show case study outputs, data tables, and literature sources to owner).

---

### Phase 5: Database & REST API
- **Goal**: Build persistent storage and asynchronous job management web API.
- **Tasks**:
  1. SQLAlchemy 2 async models: `Job`, `Protein`, `Molecule`, `Result`, `QCReport`, `ReliabilityResult`, `Validation`.
  2. Pydantic v2 schemas for all request/response objects with strict validation.
  3. FastAPI endpoints under `/api/v1`:
     - Jobs creation, protein upload/PDB fetch, ligand ingestion, pocket suggestion, box definition, settings update.
     - QC invocation and retrieval with acknowledgment gating.
     - Job execution, live status/progress, cancellation, ranked results, pose streaming, data exports.
     - Example loading endpoints (`/examples`, `/examples/1hsg/load`).
  4. Storage layer abstraction (`LocalStorageManager` under `/data` volume).
  5. Comprehensive endpoint tests with `pytest` using SQLite test database.
- **Acceptance Criteria**: 100% route coverage passing tests, strict error JSON schema, zero raw 500 tracebacks leaked.
- **Approval Gate**: No.

---

### Phase 6: Asynchronous Task Queue (Celery + Redis)
- **Goal**: Build distributed task execution for non-blocking computation.
- **Tasks**:
  1. Celery configuration with Redis broker and backend.
  2. Granular task design:
     - Molecule preparation task.
     - Individual molecule docking task.
     - Job aggregator task (updates live progress after each molecule finishes).
     - Reliability check sub-tasks for top-K candidates.
  3. Robustness controls: per-molecule timeout (10 min), failure isolation (one bad SMILES does not abort the job), cancellation handler, 7-day retention cleaner.
  4. Learning chapter: `10-backend-and-queue.md`.
- **Acceptance Criteria**: API starts job and returns immediately; progress polls smoothly; individual malformed ligands fail gracefully without killing sibling jobs.
- **Approval Gate**: No.

---

### Phase 7: Minimalist Scientific Frontend
- **Goal**: Construct the data-first, lab-notebook web interface following Section 12 design specifications.
- **Tasks**:
  1. Vite + React 18 + TypeScript setup with self-hosted Newsreader (serif), IBM Plex Mono (monospace), and plain CSS variables. Zero Tailwind, zero component libraries.
  2. Single progressive scrolling layout with sticky section index:
     - Section 1: Inputs & Biological Target Card.
     - Section 2: Search Box & Pocket Inspector.
     - Section 3: Pre-Docking QC Checklist (with explicit ACK toggle for failures).
     - Section 4: Parameters & Docking Settings.
     - Section 5: Ranked Results Table, Residue Contact Strip, and Drug-Likeness Metrics.
  3. Reliability results table with inline SVG seed score distribution dot plots.
  4. TanStack Query polling state, Zustand store for job configurations.
- **Acceptance Criteria**: Clean, typography-driven UI renders smoothly; loads example, triggers QC, and presents results within one unified page.
- **Approval Gate**: No.

---

### Phase 8: 3D Molecular Visualization (Mol*)
- **Goal**: Integrate Mol* for structural inspection of complexes, pockets, and poses.
- **Tasks**:
  1. Self-contained Mol* viewer wrapper component in hairline paper frame.
  2. Multi-representation rendering: receptor cartoon, docked ligand sticks, crystal ligand ghost pose, translucent search box wireframe, dashed interaction contacts.
  3. Interactive binding box controller with numeric inputs and dynamic updates.
  4. Learning chapter: `11-frontend-and-viewer.md`.
- **Acceptance Criteria**: Clicking result table rows updates ligand pose; displays contacts; bounding box matches backend coordinates.
- **Approval Gate**: No.

---

### Phase 9: Reports & Scientific Exports
- **Goal**: Provide reproducible data exports and publication-grade reports.
- **Tasks**:
  1. Formats: CSV summary, SDF poses archive, settings JSON (full reproduction specification), ZIP bundle.
  2. PDF generation via WeasyPrint: clean technical lab report rendering target card, QC audit log, redocking metrics, reliability panel, and software version manifest.
- **Acceptance Criteria**: Re-docking from exported settings JSON and seed reproduces numerical scores identically; PDF generates cleanly.
- **Approval Gate**: No.

---

### Phase 10: Final Polish, Documentation, & Deployment
- **Goal**: System hardening, end-to-end testing, documentation finalization, and release tagging.
- **Tasks**:
  1. End-to-end integration tests with Playwright (Load Example $\rightarrow$ QC $\rightarrow$ Run $\rightarrow$ Export).
  2. Accessibility (WCAG AA) and mobile responsive table scrolling audit.
  3. Compile complete documentation: `README.md`, `ARCHITECTURE.md`, `SCIENCE.md`, `LIMITATIONS.md`, `GLOSSARY.md`, `INTERVIEW-QA.md`.
  4. Production container build and deployment configuration.
  5. Git tag `v1.0.0`.
- **Acceptance Criteria**: Clean audit, all tests green, stranger can load example on live URL and inspect results.
- **Approval Gate**: **YES** (final owner review).

---

## 3. Key Assumptions

1. **Hardware & Operating System**: Development machine runs Windows 11. Because Docker Desktop / WSL is not installed locally on this host, we adhere strictly to **Execution Protocol Rule 4**: use the local Conda / Micromamba environment (`docklab`) derived from `environment.yml` for local execution and validation, while keeping Dockerfile and docker-compose assets fully production-ready for deployment.
2. **Computational Scope**: Vina exhaustiveness defaults to 8 for balanced web jobs, and 4 for rapid unit testing/preflight verification. All multi-seed calculations use explicit seed values ($1, 2, 3, \dots$).
3. **External Dependencies**: Direct internet access is utilized only for initial retrieval of PDB structures, PubChem SMILES, and UniProt records. Once fetched, files are cached locally in `examples/` so the system operates completely offline at runtime.
4. **Licensing & Tooling**: AutoDock Vina (Apache 2.0), RDKit (BSD), Open Babel (GPL), Meeko (LGPL), PDBFixer (MIT), ProLIF (Apache 2.0) are fully compatible with open academic/research distribution.

---

## 4. Identified Scientific & Engineering Risks and Mitigations

| Risk | Impact | Mitigation Strategy |
| :--- | :--- | :--- |
| **Vina Stochastic Search Variance** | Poses may vary between runs with different seeds or CPU thread counts. | Pin `cpu=1` and `seed` for deterministic regression tests; enforce tolerance bands ($\pm 0.3 \text{ kcal/mol}$) in regression assertions; feature this exact behavior as the scientific motivation for Feature 2 (Reliability Check). |
| **Receptor Preparation Failures with Meeko** | Some PDB structures with non-standard residues fail Meeko receptor formatting. | Implement an automated fallback pipeline: attempt Meeko first; if conversion throws an exception, log the event in `warnings` and convert via Open Babel (`obabel -ipdb ... -opdbqt -xr`). |
| **Symmetry-Aware RMSD Discrepancies** | Chemically equivalent symmetric atoms (e.g. carboxylate oxygens, phenyl carbons) yield erroneously inflated RMSD if calculated naively. | Utilize `spyrmsd` or RDKit's `GetBestRMS` which performs graph isomorphism matching to ensure true minimum symmetric RMSD. |
| **Bridging Water & Metal Charge Handling** | Vina does not explicitly model solvent molecules or divalent coordination geometry by default. | Surface explicit QC warnings when structural waters or metal ions exist in the active site; educate the user in the QC checklist and learning guide. |
| **PDB Coordinate Gaps Near Pocket** | Unmodeled disordered loops near the binding pocket alter electrostatics and steric boundaries. | QC check parses REMARK 465 and checks distance of missing residues to box center; marks as `FAIL` if within 8 Å of the pocket, requiring explicit user acknowledgement before proceeding. |

---

## 5. Learning Guide Chapter Mapping

Every chapter follows the strict 8-part template defined in Section 11 and Section 19.6:

| Chapter | Title | Phase Delivered |
| :--- | :--- | :--- |
| `00` | `00-big-picture.md` | Phase 0 |
| `01` | `01-environment-and-tools.md` | Phase 0 |
| `13` | `13-protein-structure-and-binding-basics.md` | Phase 0 |
| `14` | `14-hiv-protease-and-drug-resistance.md` | Phase 0 |
| `02` | `02-protein-preparation.md` | Phase 1 |
| `03` | `03-ligand-preparation.md` | Phase 1 |
| `04` | `04-search-box.md` | Phase 1 |
| `05` | `05-docking-with-vina.md` | Phase 1 |
| `06` | `06-redocking-and-rmsd.md` | Phase 1 |
| `15` | `15-kinases-and-egfr.md` | Phase 1 |
| `16` | `16-reading-a-docking-score-correctly.md` | Phase 1 |
| `17` | `17-drug-likeness-and-ligand-chemistry.md` | Phase 1 |
| `18` | `18-interaction-fingerprints.md` | Phase 1 |
| `07` | `07-qc-engine.md` | Phase 2 |
| `08` | `08-reliability-check.md` | Phase 3 |
| `09` | `09-roc-auc-enrichment.md` | Phase 4 |
| `12` | `12-benchmark.md` | Phase 4 |
| `19` | `19-case-study-walkthrough.md` | Phase 4.5 |
| `10` | `10-backend-and-queue.md` | Phase 6 |
| `11` | `11-frontend-and-viewer.md` | Phase 8 |
| `Glossary` | `GLOSSARY.md` | Phase 10 |
| `QA` | `INTERVIEW-QA.md` | Phase 10 |

