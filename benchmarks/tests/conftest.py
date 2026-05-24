"""Make ``cbench`` importable from the benchmarks/ root when running pytest."""

import sys
from pathlib import Path

_BENCHMARKS_ROOT = Path(__file__).resolve().parent.parent
if str(_BENCHMARKS_ROOT) not in sys.path:
    sys.path.insert(0, str(_BENCHMARKS_ROOT))
