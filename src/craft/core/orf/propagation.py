"""Reference-isoform ORF propagation.

Project the parent transcript's CDS coordinates onto the novel isoform where structure
is preserved. At divergence, predict downstream frame impact: continuation, frameshift,
premature stop, alternative stop.
"""

import pyranges as pr


def propagate(
    isoforms: pr.PyRanges,
    reference: pr.PyRanges,
    completeness: pr.PyRanges,
) -> pr.PyRanges:
    """Propagate CDS coordinates from parent reference transcripts onto novel isoforms."""
    raise NotImplementedError
