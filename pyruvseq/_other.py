"""RUVr/RUVs matrix methods; see pinned RUVSeq source in manifest.yaml.

RUVSeq authors: Davide Risso and Sandrine Dudoit; Artistic-2.0.
"""
import warnings
import numpy as np
import pandas as pd
from scipy import sparse


def _prepare(x, controls, k, epsilon, tolerance, is_log):
    labels = x if isinstance(x, pd.DataFrame) else None
    values = x.to_numpy() if labels is not None else x
    values = values.toarray() if sparse.issparse(values) else np.asarray(values)
    if (values.ndim != 2 or min(values.shape) == 0 or
            not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values) or
            not np.isfinite(values).all()):
        raise ValueError('x must be a finite genes-by-samples matrix.')
    if isinstance(k, (bool, np.bool_)) or not isinstance(k, (int, np.integer)) or not 0 < k < values.shape[1]:
        raise ValueError('k must be a positive integer less than the sample count.')
    if not np.isfinite(epsilon) or epsilon <= 0 or not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError('Require positive epsilon and nonnegative tolerance.')
    if not is_log and (np.any(values < 0) or not np.allclose(values, np.rint(values), rtol=0, atol=1e-8)):
        raise ValueError('Use raw integer counts or set is_log=True.')
    selected = np.asarray(controls)
    if selected.ndim != 1 or not len(selected):
        raise ValueError('Select a nonempty one-dimensional set of controls.')
    if labels is not None and (not labels.index.is_unique or not labels.columns.is_unique):
        raise ValueError('Gene and sample labels must be unique.')
    if selected.dtype.kind == 'b':
        if len(selected) != values.shape[0]:
            raise ValueError('Control mask length mismatch.')
        idx = np.flatnonzero(selected)
    elif labels is not None:
        idx = labels.index.get_indexer(selected)
    elif selected.dtype.kind in 'iu':
        idx = selected
    else:
        raise ValueError('Array controls must be integer indices or a boolean mask.')
    if not len(idx) or len(np.unique(idx)) != len(idx) or np.any(idx < 0) or np.any(idx >= values.shape[0]):
        raise ValueError('Controls must be distinct valid genes.')
    y = values.astype(float).T
    return (y if is_log else np.log(y + epsilon)), idx, labels


def _rank(d, k, tolerance):
    effective = min(k, int(np.count_nonzero(d > tolerance)))
    if not effective:
        raise ValueError('No nonzero control subspace remains.')
    if effective < k:
        warnings.warn(f'Rank reduces k from {k} to {effective}.', UserWarning, stacklevel=3)
    return effective


def _result(w, corrected, labels, epsilon, is_log, round_counts):
    if not is_log:
        with np.errstate(over='raise'):
            corrected = np.exp(corrected) - epsilon
        if round_counts:
            corrected = np.maximum(np.rint(corrected), 0)
    corrected = corrected.T
    if labels is not None:
        corrected = pd.DataFrame(corrected, index=labels.index, columns=labels.columns)
        w = pd.DataFrame(w, index=labels.columns, columns=[f'W_{i+1}' for i in range(w.shape[1])])
    return {'W': w, 'normalized_counts': corrected, 'k_used': w.shape[1]}


def ruvr(x, controls, k, residuals, *, center=True, round_counts=True,
         epsilon=1., tolerance=1e-8, is_log=False):
    """Estimate unwanted factors from caller-supplied genes-by-samples residuals.

    Residuals should come from an appropriate first-pass biological model; this
    function does not fit that model. Use returned W with original counts in DE.
    """
    y, idx, labels = _prepare(x, controls, k, epsilon, tolerance, is_log)
    if isinstance(residuals, pd.DataFrame) and labels is not None:
        if not residuals.index.equals(labels.index) or not residuals.columns.equals(labels.columns):
            raise ValueError('Residual labels must match expression axes in order.')
    e = np.asarray(residuals, dtype=float).T
    if e.shape != y.shape or not np.isfinite(e).all():
        raise ValueError('Residuals must be finite and have the same shape as x.')
    if center:
        e = e - e.mean(axis=0)
    u, d, _ = np.linalg.svd(e[:, idx], full_matrices=False)
    w = u[:, :_rank(d, k, tolerance)]
    return _result(w, y - w @ (w.T @ y), labels, epsilon, is_log, round_counts)


def make_groups(groups):
    """Return zero-based replicate indices, padded with -1 (R uses one-based)."""
    groups = np.asarray(groups)
    if groups.ndim != 1 or not len(groups) or pd.isna(groups).any():
        raise ValueError('groups must contain one nonmissing label per sample.')
    rows = [np.flatnonzero(groups == value) for value in sorted(pd.unique(groups))]
    result = np.full((len(rows), max(map(len, rows))), -1, dtype=int)
    for i, row in enumerate(rows):
        result[i, :len(row)] = row
    return result


def ruvs(x, controls, k, replicate_indices, *, round_counts=True,
         epsilon=1., tolerance=1e-8, is_log=False):
    """Estimate factors using replicate groups with constant biological signal.

    replicate_indices is a matrix of zero-based sample positions padded with
    -1, as returned by make_groups. Singleton groups do not estimate variation.
    """
    y, idx, labels = _prepare(x, controls, k, epsilon, tolerance, is_log)
    groups = np.asarray(replicate_indices)
    if groups.ndim != 2 or groups.dtype.kind not in 'iu' or np.any(groups < -1) or np.any(groups >= len(y)):
        raise ValueError('replicate_indices must contain valid zero-based positions or -1.')
    centered = []
    for row in groups:
        row = row[row >= 0]
        if len(np.unique(row)) != len(row):
            raise ValueError('Replicate rows cannot repeat a sample.')
        if len(row) >= 2:
            centered.extend(y[row] - y[row].mean(axis=0))
    if not centered:
        raise ValueError('At least one replicate group needs two samples.')
    yc = np.asarray(centered)
    yc = yc[yc.sum(axis=1) != 0]  # preserve upstream row filtering semantics
    if not len(yc):
        raise ValueError('Replicates provide no variation.')
    _, d, vt = np.linalg.svd(yc, full_matrices=False)
    effective = _rank(d, k, tolerance)
    alpha = d[:effective, None] * vt[:effective]
    ac = alpha[:, idx]
    if np.linalg.matrix_rank(ac) < effective:
        raise ValueError('Control loadings cannot identify all requested factors.')
    w = y[:, idx] @ np.linalg.solve(ac @ ac.T, ac).T
    return _result(w, y - w @ alpha, labels, epsilon, is_log, round_counts)
