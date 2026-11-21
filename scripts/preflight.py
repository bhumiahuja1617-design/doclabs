#!/usr/bin/env python3
"""
scripts/preflight.py
DockLab Preflight Verification Script.

Inspects environment, imports all critical computational chemistry and web
packages, invokes core external binaries (vina, fpocket, obabel), outputs
their exact versions, and handles platform differences honestly.
"""

import sys
import os
import subprocess
import importlib
from importlib import metadata
from pathlib import Path

REQUIRED_PYTHON_PACKAGES = [
    # Chemistry & Structural Modeling
    ("rdkit", "rdkit.__version__"),
    ("openbabel", "openbabel.__version__"),
    ("meeko", "meeko.__version__"),
    ("pdbfixer", "pdbfixer.__version__"),
    ("openmm", "openmm.__version__"),
    ("Bio", "Bio.__version__"),  # Biopython
    ("prolif", "prolif.__version__"),
    ("spyrmsd", "spyrmsd.__version__"),
    ("vina", "vina.__version__"),
    ("weasyprint", "weasyprint.__version__"),
    ("pydantic_settings", "pydantic_settings.__version__"),
    ("MDAnalysis", "MDAnalysis.__version__"),
    
    # Scientific & ML
    ("numpy", "numpy.__version__"),
    ("scipy", "scipy.__version__"),
    ("sklearn", "sklearn.__version__"),

    # Web & Application Infrastructure
    ("fastapi", "fastapi.__version__"),
    ("uvicorn", "uvicorn.__version__"),
    ("celery", "celery.__version__"),
    ("redis", "redis.__version__"),
    ("sqlalchemy", "sqlalchemy.__version__"),
    ("alembic", "alembic.__version__"),
    ("pydantic", "pydantic.__version__"),
]

def check_package(module_name: str, version_expr: str) -> tuple[bool, str]:
    try:
        mod = importlib.import_module(module_name)
        try:
            ver = eval(version_expr, {module_name: mod})
        except Exception:
            distributions = metadata.packages_distributions().get(module_name, [])
            ver = next((metadata.version(name) for name in distributions), getattr(mod, "__version__", "unknown"))
        return True, str(ver)
    except Exception as exc:
        return False, str(exc)

def check_binary(bin_name: str, default_cmd: list[str]) -> tuple[str, str]:
    # Check candidates: ./bin, conda env Library/bin, conda env Scripts, and system PATH
    root_dir = Path(__file__).resolve().parent.parent
    local_bin = root_dir / "bin"
    python_dir = Path(sys.executable).parent
    conda_library_bin = python_dir / "Library" / "bin"
    conda_scripts = python_dir / "Scripts"

    search_dirs = [local_bin, conda_library_bin, conda_scripts]
    
    candidate_cmds = []
    for d in search_dirs:
        exe_path = d / f"{bin_name}.exe" if sys.platform == "win32" else d / bin_name
        if exe_path.exists():
            candidate_cmds.append([str(exe_path)] + default_cmd[1:])
    
    candidate_cmds.append(default_cmd)

    for cmd in candidate_cmds:
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=10,
            )
            combined = (result.stdout + "\n" + result.stderr).strip()
            if result.returncode == 0 or (bin_name == "fpocket" and "fpocket" in combined.lower()):
                output_line = combined.splitlines()[0] if combined else "Available"
                return "PASS", output_line
        except FileNotFoundError:
            continue
        except Exception as exc:
            return "FAIL", str(exc)

    if bin_name == "fpocket" and sys.platform == "win32":
        return "WARN", "UNAVAILABLE on Windows (Linux/Docker only). Pocket detection degrades gracefully."
    
    return "FAIL", "Executable not found on system PATH or ./bin"

def main() -> int:
    print("=" * 70)
    print("  DockLab Phase 0 Preflight Verification")
    print(f"  Platform: {sys.platform} | Python: {sys.version.split()[0]}")
    print(f"  Executable: {sys.executable}")
    print("=" * 70)
    
    critical_failures = 0

    print("\n--- [1] Checking Required Python Packages ---")
    for mod_name, ver_expr in REQUIRED_PYTHON_PACKAGES:
        ok, detail = check_package(mod_name, ver_expr)
        status = "PASS" if ok else "FAIL"
        print(f"[{status:4s}] {mod_name:<16}: {detail}")
        if not ok:
            critical_failures += 1

    print("\n--- [2] Checking Required Scientific Binaries ---")
    binaries = [
        ("vina", ["vina", "--help"]),
        ("fpocket", ["fpocket", "-h"]),
        ("obabel", ["obabel", "-V"]),
    ]
    for bin_name, cmd in binaries:
        status, detail = check_binary(bin_name, cmd)
        print(f"[{status:4s}] {bin_name:<16}: {detail}")
        if status == "FAIL":
            critical_failures += 1

    print("\n" + "=" * 70)
    if critical_failures == 0:
        print(">> ALL CRITICAL PREFLIGHT CHECKS PASSED. Ready for science core. <<")
        print("=" * 70)
        return 0
    else:
        print(f">> PREFLIGHT FAILED: {critical_failures} critical check(s) did not pass. <<")
        print("=" * 70)
        return 1

if __name__ == "__main__":
    sys.exit(main())
