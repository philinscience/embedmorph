"""Procrustes alignment between two point clouds representing the same cells."""

from __future__ import annotations

import numpy as np


def procrustes_align(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotate, reflect, and scale ``source`` onto ``target``'s frame.

    ``source`` and ``target`` must be ``(n_cells, 2)`` arrays whose rows are the
    same cells in the same order (e.g. two ``.obsm`` embeddings of one AnnData).
    Returns ``source`` transformed by the similarity transform (translation +
    rotation/reflection + isotropic scale) that minimizes summed squared
    distance to ``target``, so a straight-line interpolation between the
    result and ``target`` moves cells coherently instead of across the whole
    canvas.
    """
    if source.shape != target.shape:
        raise ValueError(f"shape mismatch: source {source.shape} vs target {target.shape}")
    if source.shape[1] != 2:
        raise ValueError(f"expected 2D coordinates, got shape {source.shape}")

    source_mean = source.mean(axis=0)
    target_mean = target.mean(axis=0)
    source_c = source - source_mean
    target_c = target - target_mean

    u, s, vt = np.linalg.svd(source_c.T @ target_c)
    rotation = u @ vt

    source_norm_sq = float((source_c**2).sum())
    scale = s.sum() / source_norm_sq if source_norm_sq > 0 else 1.0

    return scale * (source_c @ rotation) + target_mean
