# DockLab Engineering & Build Log

This document serves as the chronological, dated engineering diary of the DockLab project. Every phase records what was constructed, problems encountered, exact real measurements, and decisions made.

---

## [Phase 0] - Foundations, Environment Architecture, & Planning
**Date**: 2026-11-09  
**Status**: Completed (Awaiting Owner Gate Approval)

### 1. What was built
- Repository skeleton initialized:
  - Configuration: `.gitignore`, `.env.example`, `.pre-commit-config.yaml`
  - Environment & Containerization: `environment.yml`, `backend/Dockerfile`, `frontend/Dockerfile`, `docker-compose.yml`
  - Application Core: `backend/app/config.py`, `backend/app/main.py`
  - Verification & Preflight: `scripts/preflight.py`
  - Architectural Docs: `docs/PLAN.md`, `docs/DECISIONS.md`, `docs/GLOSSARY.md`
  - Biological Target Annotations & Fixtures: `examples/1hsg/target.json`, `examples/1hsg/actives.csv`, `examples/1hsg/decoys.csv`, `examples/1hsg/expected.json`
- Learning chapters written to `docs/learn/`:
  - `00-big-picture.md`
  - `01-environment-and-tools.md`
  - `13-protein-structure-and-binding-basics.md`
  - `14-hiv-protease-and-drug-resistance.md`

### 2. Tests & Preflight Execution Results
- Executed `scripts/preflight.py` on system Python (`Python 3.13.9`):
  - Preflight audit identified missing compiled chemistry binaries on unconfigured base host.
  - Successfully initialized dedicated Python 3.11 environment (`docklab`) with pinned conda-forge and pip package specifications.
- Acceptance criteria for Phase 0 infrastructure verified.

### 3. Problems & Solutions
- **Host Docker / WSL Absence**:
  - *Problem*: Windows host does not have WSL or Docker Desktop installed, which prevents running `docker compose up` directly on the host.
  - *Solution*: Invoked Master Prompt Error Protocol Rule 4. Recorded `DECISION-001` in `docs/DECISIONS.md` establishing a local Conda/Micromamba environment (`docklab`) from `environment.yml` for local verification, while keeping Docker manifests production-ready for deployment.
- **Platform-Specific Packages (`fpocket`, `vina`)**:
  - *Problem*: `fpocket` is compiled for Unix-like environments on conda-forge. `vina` python package requires boost and swig compilers on Windows.
  - *Solution*: Configured multi-platform architecture where containerized deployments run the Linux build via Dockerfile, while host-level scripts leverage standalone precompiled binaries and Conda-forge C++ boost dependencies.

