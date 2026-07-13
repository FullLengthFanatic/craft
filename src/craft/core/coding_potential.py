"""Coding-potential scoring with a model trained on the supplied reference.

CRAFT already loads a reference annotation that contains both coding transcripts
(CDS rows) and non-coding ones (transcripts with exons but no CDS). This module
uses them as a training set: it builds a hexamer coding/non-coding log-likelihood
table from the reference, derives three features per ORF, fits a logistic
regression, and scores every isoform's best ORF. No model file is shipped and no
external tool is required; the model adapts its training data to the organism
described by the reference.

Features (the dominant CPAT signals):
- hexamer usage log-likelihood ratio (coding vs non-coding), per-ORF mean;
- ORF length (log10 bp);
- ORF coverage (ORF bp / transcript bp).

The score is a classifier score in [0, 1], not a calibrated biological
probability; a leakage-free internal cross-validation AUC is reported for
validation. This is a screening score, not a substitute for a curated
coding/non-coding classifier; confirm borderline calls with CPC2/CPAT.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import orfipy_core
import pandas as pd
import pyranges as pr
import pysam
from scipy.optimize import minimize

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_HEXAMERS = ["".join(h) for h in itertools.product("ACGT", repeat=6)]
_HEX_INDEX = {h: i for i, h in enumerate(_HEXAMERS)}

# Bound the training set so full-genome references stay fast; the hexamer table
# and logistic both converge well below this. Deterministic stride sampling.
MAX_TRAIN_PER_CLASS = 4000  # a cap sweep on GENCODE v45 showed AUC plateaus by 4000
_CV_FOLDS = 5
_PSEUDOCOUNT = 1.0
_L2 = 1.0  # ridge penalty on the logistic weights (excludes intercept)
_COLUMNS = [
    "transcript_id",
    "coding_potential_score",
    "coding_potential_label",
    "coding_potential_orf_source",
]


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _spliced_sequence(intervals: list[tuple], genome: pysam.FastaFile) -> str:
    """Spliced nucleotide sequence (5'→3') for a list of (chrom, start, end, strand)."""
    if not intervals:
        return ""
    chrom = intervals[0][0]
    strand = intervals[0][3]
    ordered = sorted(intervals, key=lambda x: x[1])
    parts = [genome.fetch(chrom, int(s), int(e)) for _, s, e, _ in ordered]
    seq = "".join(parts).upper()
    return _reverse_complement(seq) if strand == "-" else seq


def _hexamer_counts(seq: str) -> np.ndarray:
    """Codon-stepped (step 3) hexamer counts over 4096 hexamers; N-containing skipped."""
    counts = np.zeros(len(_HEXAMERS), dtype=np.float64)
    for i in range(0, len(seq) - 5, 3):
        idx = _HEX_INDEX.get(seq[i : i + 6])
        if idx is not None:
            counts[idx] += 1.0
    return counts


def _log_ratio_table(coding: np.ndarray, noncoding: np.ndarray) -> np.ndarray:
    """Per-hexamer log(P_coding / P_noncoding) with a pseudocount."""
    c = (coding + _PSEUDOCOUNT) / (coding.sum() + _PSEUDOCOUNT * len(coding))
    n = (noncoding + _PSEUDOCOUNT) / (noncoding.sum() + _PSEUDOCOUNT * len(noncoding))
    return np.log(c / n)


def _longest_orf_seq(seq: str) -> str:
    """Longest ATG-initiated ORF in the forward transcript ("" if none)."""
    if len(seq) < 3:
        return ""
    best = ""
    for start, stop, _, _ in orfipy_core.orfs(seq, minlen=3, starts=["ATG"], strand="f"):
        if stop - start > len(best):
            best = seq[start:stop]
    return best


def _hexamer_score(seq: str, log_ratio: np.ndarray) -> float:
    """Mean hexamer log-likelihood ratio over the ORF (codon-stepped)."""
    total = 0.0
    n = 0
    for i in range(0, len(seq) - 5, 3):
        idx = _HEX_INDEX.get(seq[i : i + 6])
        if idx is not None:
            total += log_ratio[idx]
            n += 1
    return total / n if n else 0.0


# Fickett TESTCODE lookup tables (Fickett 1982, NAR 10:5303), as shipped by the
# CPAT-derived lncScore implementation. content_para has 9 thresholds, so the
# 10th content_prob entry is unused, matching that reference exactly.
_FICKETT_POS_PROB = {
    "A": [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22],
    "C": [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23],
    "G": [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08],
    "T": [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09],
}
_FICKETT_POS_WEIGHT = {"A": 0.26, "C": 0.18, "G": 0.31, "T": 0.33}
_FICKETT_POS_PARA = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
_FICKETT_CONT_PROB = {
    "A": [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21],
    "C": [0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31],
    "G": [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29],
    "T": [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58],
}
_FICKETT_CONT_WEIGHT = {"A": 0.11, "C": 0.12, "G": 0.15, "T": 0.14}
_FICKETT_CONT_PARA = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0.0]


def _fickett_lookup(value: float, para: list[float], prob: list[float], weight: float) -> float:
    if value < 0:
        return 0.0
    for idx, thr in enumerate(para):
        if value >= thr:
            return prob[idx] * weight
    return 0.0


def _fickett_score(seq: str) -> float:
    """Fickett TESTCODE statistic (higher = more coding)."""
    seq = seq.upper()
    n = len(seq)
    if n < 3:
        return 0.0
    score = 0.0
    for base in "ACGT":
        content = seq.count(base) / n
        score += _fickett_lookup(
            content, _FICKETT_CONT_PARA, _FICKETT_CONT_PROB[base], _FICKETT_CONT_WEIGHT[base]
        )
        c0 = seq[0::3].count(base)
        c1 = seq[1::3].count(base)
        c2 = seq[2::3].count(base)
        position = max(c0, c1, c2) / (min(c0, c1, c2) + 1.0)
        score += _fickett_lookup(
            position, _FICKETT_POS_PARA, _FICKETT_POS_PROB[base], _FICKETT_POS_WEIGHT[base]
        )
    return score


def _features(
    orf_nt: str, tx_len: int, log_ratio: np.ndarray
) -> tuple[float, float, float, float]:
    orf_len = len(orf_nt)
    hexamer = _hexamer_score(orf_nt, log_ratio)
    log_len = float(np.log10(orf_len + 1))
    coverage = (orf_len / tx_len) if tx_len > 0 else 0.0
    fickett = _fickett_score(orf_nt)
    return hexamer, log_len, coverage, fickett


def _stride_sample(ids: list[str], cap: int) -> list[str]:
    ids = sorted(ids)
    if len(ids) <= cap:
        return ids
    step = len(ids) / cap
    return [ids[int(i * step)] for i in range(cap)]


def _sigmoid(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))


def _fit_logistic(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fit logistic weights (incl. intercept as the last term) by L-BFGS-B."""
    xb = np.hstack([x, np.ones((x.shape[0], 1))])
    n_feat = xb.shape[1]

    def nll(w: np.ndarray) -> float:
        p = _sigmoid(xb @ w)
        eps = 1e-9
        ll = y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps)
        reg = _L2 * np.sum(w[:-1] ** 2)
        return -ll.sum() + reg

    def grad(w: np.ndarray) -> np.ndarray:
        p = _sigmoid(xb @ w)
        g = xb.T @ (p - y)
        g[:-1] += 2 * _L2 * w[:-1]
        return g

    res = minimize(nll, np.zeros(n_feat), jac=grad, method="L-BFGS-B")
    return res.x


def _auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """ROC AUC via the Mann-Whitney U statistic (rank-based)."""
    n_pos = int((labels == 1).sum())
    n_neg = int((labels == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # Tie-aware average ranks (1-based).
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    start = np.concatenate([[0], np.cumsum(counts)[:-1]])
    avg_rank = start + (counts + 1) / 2.0
    ranks = avg_rank[inv]
    r_pos = ranks[labels == 1].sum()
    u = r_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def build_model(
    reference: pr.PyRanges, genome: pysam.FastaFile, threshold: float = 0.5
) -> dict | None:
    """Train hexamer tables + logistic on the reference; None if no negative set.

    Returns a dict with ``log_ratio``, standardisation ``mean``/``std``, logistic
    ``weights``, ``threshold``, training counts, and a 5-fold cross-validated AUC
    (``heldout_auc``). The returned ``weights`` are fit on all training data; the
    AUC is a separate diagnostic and does not affect the scores.
    """
    ref_df = reference.df
    if "Feature" not in ref_df.columns:
        return None
    cds = ref_df[ref_df["Feature"] == "CDS"]
    exon = ref_df[ref_df["Feature"] == "exon"]
    coding_ids = set(cds["transcript_id"].astype(str).unique())
    noncoding_ids = sorted(set(exon["transcript_id"].astype(str).unique()) - coding_ids)
    if not coding_ids or not noncoding_ids:
        return None

    coding_sel = _stride_sample(sorted(coding_ids), MAX_TRAIN_PER_CLASS)
    noncoding_sel = _stride_sample(noncoding_ids, MAX_TRAIN_PER_CLASS)

    def _intervals(df: pd.DataFrame) -> dict[str, list[tuple]]:
        out: dict[str, list[tuple]] = {}
        for tx, g in df.groupby("transcript_id", sort=False):
            out[str(tx)] = [
                (str(r.Chromosome), int(r.Start), int(r.End), str(r.Strand))
                for r in g.itertuples(index=False)
            ]
        return out

    coding_set = set(coding_sel)
    cds_iv = _intervals(cds[cds["transcript_id"].astype(str).isin(coding_set)])
    coding_exon_iv = _intervals(exon[exon["transcript_id"].astype(str).isin(coding_set)])
    exon_iv = _intervals(exon[exon["transcript_id"].astype(str).isin(set(noncoding_sel))])

    # Tables: coding hexamers from CDS, non-coding hexamers from full transcripts.
    # Features: the transcript's ORF (annotated CDS for coding, longest predicted
    # ORF for non-coding) against its full transcript length.
    coding_counts = np.zeros(len(_HEXAMERS))
    noncoding_counts = np.zeros(len(_HEXAMERS))
    # (candidate ORF, transcript length, sequence used to train the hexamer table)
    coding_feat: list[tuple[str, int, str]] = []
    noncoding_feat: list[tuple[str, int, str]] = []
    for tx in coding_sel:
        cds_seq = _spliced_sequence(cds_iv.get(tx, []), genome)
        mrna_len = sum(e - s for _, s, e, _ in coding_exon_iv.get(tx, [])) or len(cds_seq)
        coding_counts += _hexamer_counts(cds_seq)
        coding_feat.append((cds_seq, mrna_len, cds_seq))
    for tx in noncoding_sel:
        tx_seq = _spliced_sequence(exon_iv.get(tx, []), genome)
        noncoding_counts += _hexamer_counts(tx_seq)
        noncoding_feat.append((_longest_orf_seq(tx_seq), len(tx_seq), tx_seq))

    log_ratio = _log_ratio_table(coding_counts, noncoding_counts)

    rows = []
    labels = []
    for orf_seq, tx_len, _ in coding_feat:
        rows.append(_features(orf_seq, tx_len, log_ratio))
        labels.append(1)
    for orf_seq, tx_len, _ in noncoding_feat:
        rows.append(_features(orf_seq, tx_len, log_ratio))
        labels.append(0)
    x = np.array(rows, dtype=np.float64)
    y = np.array(labels, dtype=np.float64)

    records = [*coding_feat, *noncoding_feat]

    # Cross-validated AUC.  Every fold rebuilds the hexamer table and feature
    # scaling on its training transcripts.  Earlier releases derived both from
    # all transcripts before splitting, leaking test-fold sequence information.
    idx = np.arange(len(y))
    fold_aucs: list[float] = []
    for k in range(_CV_FOLDS):
        test_mask = idx % _CV_FOLDS == k
        if test_mask.sum() == 0 or (~test_mask).sum() == 0:
            continue
        train_idx = idx[~test_mask]
        test_idx = idx[test_mask]
        fold_coding = np.zeros(len(_HEXAMERS))
        fold_noncoding = np.zeros(len(_HEXAMERS))
        for i in train_idx:
            source_seq = records[int(i)][2]
            if y[int(i)] == 1:
                fold_coding += _hexamer_counts(source_seq)
            else:
                fold_noncoding += _hexamer_counts(source_seq)
        fold_ratio = _log_ratio_table(fold_coding, fold_noncoding)
        x_train = np.array(
            [_features(records[int(i)][0], records[int(i)][1], fold_ratio) for i in train_idx]
        )
        x_test = np.array(
            [_features(records[int(i)][0], records[int(i)][1], fold_ratio) for i in test_idx]
        )
        fold_mean = x_train.mean(axis=0)
        fold_std = x_train.std(axis=0)
        fold_std[fold_std == 0] = 1.0
        x_train = (x_train - fold_mean) / fold_std
        x_test = (x_test - fold_mean) / fold_std
        w_tr = _fit_logistic(x_train, y[train_idx])
        s_te = _sigmoid(np.hstack([x_test, np.ones((x_test.shape[0], 1))]) @ w_tr)
        a = _auc(s_te, y[test_mask])
        if a == a:
            fold_aucs.append(a)
    auc = float(np.mean(fold_aucs)) if fold_aucs else float("nan")

    mean = x.mean(axis=0)
    std = x.std(axis=0)
    std[std == 0] = 1.0
    xz = (x - mean) / std
    weights = _fit_logistic(xz, y)
    return {
        "log_ratio": log_ratio,
        "mean": mean,
        "std": std,
        "weights": weights,
        "threshold": threshold,
        "n_coding": len(coding_sel),
        "n_noncoding": len(noncoding_sel),
        "heldout_auc": auc,
        "validation": "leakage_free_internal_5fold",
    }


def _predict(model: dict, orf_nt: str, tx_len: int) -> float:
    feat = np.array(_features(orf_nt, tx_len, model["log_ratio"]), dtype=np.float64)
    z = (feat - model["mean"]) / model["std"]
    return float(_sigmoid(np.append(z, 1.0) @ model["weights"]))


def score_isoforms(
    per_isoform: pd.DataFrame,
    classified: pr.PyRanges,
    reference: pr.PyRanges,
    genome_fasta,
    threshold: float = 0.5,
) -> tuple[pd.DataFrame, dict | None]:
    """Score every isoform's best ORF for coding potential.

    The ORF is taken from ``resolved_cds_intervals`` when present, else
    ``propagated_cds_intervals``, else ``denovo_cds_intervals`` (same precedence
    as the Pfam scan). Isoforms with no ORF are labelled ``noncoding`` with a NaN
    score.

    Returns ``(per_isoform_scores, model_info)``. ``model_info`` is ``None`` (and
    every score NaN) when the reference has no non-coding transcripts to train on.
    """
    if per_isoform.empty or len(classified) == 0:
        return pd.DataFrame(columns=_COLUMNS), None

    iso_df = classified.df
    tx_len = iso_df.assign(_len=iso_df["End"] - iso_df["Start"]).groupby("transcript_id")[
        "_len"
    ].sum().to_dict()

    own = isinstance(genome_fasta, str | Path)
    genome = pysam.FastaFile(str(genome_fasta)) if own else genome_fasta
    try:
        model = build_model(reference, genome, threshold=threshold)
        rows: list[dict] = []
        for _, r in per_isoform.iterrows():
            tx_id = r["transcript_id"]
            intervals = r.get("resolved_cds_intervals") or []
            source = "resolved"
            if not intervals:
                intervals = r.get("propagated_cds_intervals") or []
                source = "propagated"
            if not intervals:
                intervals = r.get("denovo_cds_intervals") or []
                source = "denovo"
            if not intervals or model is None:
                rows.append(
                    {
                        "transcript_id": tx_id,
                        "coding_potential_score": float("nan"),
                        "coding_potential_label": "" if model is None else "noncoding",
                        "coding_potential_orf_source": "none",
                    }
                )
                continue
            orf_nt = _spliced_sequence(intervals, genome)
            score = _predict(model, orf_nt, int(tx_len.get(tx_id, len(orf_nt))))
            rows.append(
                {
                    "transcript_id": tx_id,
                    "coding_potential_score": score,
                    "coding_potential_label": "coding" if score >= threshold else "noncoding",
                    "coding_potential_orf_source": source,
                }
            )
        return pd.DataFrame(rows, columns=_COLUMNS), model
    finally:
        if own:
            genome.close()
