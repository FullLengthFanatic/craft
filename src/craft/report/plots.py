"""Plotly figure builders for the HTML report."""

import pandas as pd
import plotly.graph_objects as go


def isoform_track(per_isoform: pd.DataFrame, gene_id: str) -> go.Figure:
    """Render an exon/CDS/domain track view for all isoforms of a single gene."""
    raise NotImplementedError
