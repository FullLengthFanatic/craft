"""Top-level HTML report assembly."""

from pathlib import Path

import pandas as pd


def render(per_isoform: pd.DataFrame, output: Path) -> None:
    """Render the interactive HTML report (summary, per-gene views, sortable table)."""
    raise NotImplementedError
