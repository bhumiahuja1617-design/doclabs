# DockLab Architecture & Scientific Decisions Record (ADR)

---

### DECISION-001: Local Environment Fallback for Windows Host
- **Date**: 2026-11-09
- **Status**: Accepted
- **Context**: The host machine is Windows 11 without Windows Subsystem for Linux (WSL) or Docker Desktop installed. Standard `docker compose up` cannot execute directly on the bare host without Docker daemon support.
- **Decision**: In strict accordance with the Master Prompt Error Protocol (Rule 4: *"If the owner's machine cannot run Docker, fall back to a local micromamba environment from environment.yml, tell the owner, and record it"*), we maintain `backend/Dockerfile`, `frontend/Dockerfile`, and `docker-compose.yml` for containerized environments, while establishing a local Conda/Micromamba environment (`docklab`) with Python 3.11 for host-level development, testing, and script execution.
- **Consequences**: Enables execution of all science CLI tools, unit tests, and preflight scripts directly on the host machine while preserving container portability for cloud deployment.

---

### DECISION-002: Python 3.11 and Pinned Conda-Forge Science Stack
- **Date**: 2026-11-09
- **Status**: Accepted
- **Context**: Structural biology tools (OpenMM, PDBFixer, Open Babel, AutoDock Vina) have complex compiled C/C++ dependencies and shared libraries that frequently break when installed via standard Windows wheel packages.
- **Decision**: Pin Python 3.11 and install primary scientific binaries/bindings via `conda-forge`. Python application and API packages (FastAPI, Celery, ProLIF, WeasyPrint) are installed via pip.
- **Consequences**: Guarantees binary compatibility across compiled chemistry extensions and OpenMM GPU/CPU platforms.

---

### DECISION-003: Meeko Receptor Preparation with Automated Open Babel Fallback
- **Date**: 2026-11-09
- **Status**: Accepted
- **Context**: Meeko (`mk_prepare_receptor.py` / Python API) is the modern AutoDock tool for PDBQT conversion, but can fail on non-standard residues or atypical crystal coordinates.
- **Decision**: Attempt receptor conversion via Meeko; if an unresolvable parsing or atom typing error occurs, automatically fall back to Open Babel (`obabel -ipdb ... -opdbqt -xr`), logging a warning in the QC report.
- **Consequences**: Zero catastrophic failures on unstandardized experimental PDB files.

---

### DECISION-004: Rigorous Separation Between Computational Repeatability and Biological Efficacy
- **Date**: 2026-11-09
- **Status**: Accepted
- **Context**: Novices frequently mistake AutoDock Vina score stability or negative scoring for true thermodynamic $K_d$ or biological inhibition.
- **Decision**: All documentation, UI banners, and report footers explicitly state: *"Docking scores and reliability labels reflect computational stability and spatial interaction estimates. They do not predict real-world binding probability or biological efficacy."*
- **Consequences**: Scientific integrity is preserved, preventing misinterpretation of docking results.


---

### DECISION-005: fpocket Platform Availability & Graceful Degradation
- **Date**: 2026-11-09
- **Status**: Accepted
- **Context**: `fpocket` is an open-source C package compiled exclusively for POSIX Unix environments (Linux and macOS). It is not distributed on Conda-forge for Windows (`win-64`) and cannot run natively on Windows without WSL or Linux containers.
- **Decision**: On native Windows environments, `fpocket` is marked as `UNAVAILABLE (Linux/Docker only)`. DockLab's pocket detection service (`box_service.py`) degrades gracefully: pocket detection via fpocket is an optional feature. Primary bounding box generation methods—specifically **crystal-ligand centroid + padding** and **user-selected active site residues / coordinates**—remain fully operational and primary. In Linux/Docker deployments, fpocket is fully functional.
- **Consequences**: No simulated or faked pocket checks; zero blocking failures on Windows host systems while retaining full pocket capability inside Docker containers.
