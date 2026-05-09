"""Make the project's source folders importable as flat modules during tests."""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

# The project root lets us write `from scripts.X import Y`.
# Other dirs are added so each script can be imported by its own bare module name
# (these scripts are designed to be run directly, not packaged).
_EXTRA_PATHS = [
    PROJECT_ROOT,
    PROJECT_ROOT / "software",
    PROJECT_ROOT / "tutorials" / "protein_design" / "improve_protein_protein_affinity" / "scripts",
]

for _p in _EXTRA_PATHS:
    if _p.exists() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
