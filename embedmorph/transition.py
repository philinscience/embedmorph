"""Morph a spatial transcriptomics scatter plot into its UMAP embedding (or any AnnData .obsm-to-.obsm pair) as an mp4."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np

from .align import procrustes_align
from .render import render_frame, reorder_for_highlight, resolve_colors

try:
    from tqdm.auto import tqdm
except ImportError:  # tqdm is a convenience, not a hard requirement

    def tqdm(iterable, **_kwargs):
        return iterable


def _ease_in_out(t: np.ndarray) -> np.ndarray:
    """Smoothstep: zero velocity at both ends, avoids abrupt starts/stops."""
    return t * t * (3 - 2 * t)


def _ease_smoother(t: np.ndarray) -> np.ndarray:
    """Smootherstep: zero *and* zero acceleration at both ends.

    ``_ease_in_out`` still has a jump in acceleration right where the hold
    ends and the morph begins (velocity is 0 there either way, but the hold
    is at a dead stop while the morph's curve is already accelerating away
    from it) -- this is what reads as an abrupt onset even though velocity
    itself is continuous. Matching acceleration too removes that jerk.
    """
    return t**3 * (t * (t * 6 - 15) + 10)


def _linear(t: np.ndarray) -> np.ndarray:
    return t


_EASINGS = {"linear": _linear, "ease_in_out": _ease_in_out, "ease_smoother": _ease_smoother}


def _view_bounds(coords: np.ndarray, pad_frac: float) -> tuple[tuple[float, float], tuple[float, float]]:
    pad_x = pad_frac * np.ptp(coords[:, 0]) or 1.0
    pad_y = pad_frac * np.ptp(coords[:, 1]) or 1.0
    xlim = (coords[:, 0].min() - pad_x, coords[:, 0].max() + pad_x)
    ylim = (coords[:, 1].min() - pad_y, coords[:, 1].max() + pad_y)
    return xlim, ylim


def _enclosing(*bounds: tuple[tuple[float, float], tuple[float, float]]):
    xmin = min(b[0][0] for b in bounds)
    xmax = max(b[0][1] for b in bounds)
    ymin = min(b[1][0] for b in bounds)
    ymax = max(b[1][1] for b in bounds)
    return (xmin, xmax), (ymin, ymax)


def _arc_controls(a: np.ndarray, b: np.ndarray, arc_strength: float) -> np.ndarray:
    """Bezier control points bulging perpendicular to each cell's own A->B segment.

    Every cell curves to the same side of its direction of travel (a 90°
    rotation of the displacement vector), so the whole cloud reads as one
    coherent flow rather than cells arcing every which way.
    """
    d = b - a
    dist = np.linalg.norm(d, axis=1, keepdims=True)
    perp = np.stack([-d[:, 1], d[:, 0]], axis=1)
    safe_dist = np.where(dist > 1e-12, dist, 1.0)
    perp_unit = perp / safe_dist
    return (a + b) / 2 + perp_unit * dist * arc_strength


def _stagger_starts(a: np.ndarray, b: np.ndarray, stagger: float) -> np.ndarray:
    """Per-cell launch offset in [0, stagger]: cells with farther to travel leave sooner.

    Everyone still reaches t=1 by the last frame (see ``_positions_at``),
    so this only staggers *when* each cell starts moving, not when it lands
    -- reads as an organic flock departing and converging, rather than every
    cell starting and stopping in perfect lockstep.
    """
    if stagger <= 0:
        return np.zeros(len(a))
    dist = np.linalg.norm(b - a, axis=1)
    max_dist = dist.max()
    norm_dist = dist / max_dist if max_dist > 0 else np.zeros_like(dist)
    return stagger * (1 - norm_dist)


def make_transition_video(
    adata,
    obsm_from: str,
    obsm_to: str,
    color: str,
    out_path: str | Path,
    *,
    n_frames: int = 130,
    fps: int = 30,
    hold_frames: int = 40,
    align: bool = True,
    dims: tuple[int, int] = (0, 1),
    point_size: float = 8.0,
    alpha: float = 0.85,
    figsize: tuple[float, float] = (8.0, 7.0),
    dpi: int = 150,
    facecolor: str = "white",
    palette: dict | str | Sequence[str] | None = None,
    show_legend: bool = True,
    easing: str = "ease_smoother",
    labels: tuple[str, str] | None = None,
    show_progress: bool = True,
    zoom_to_fit: bool = True,
    view_padding: float = 0.08,
    arc_strength: float = 0.0,
    stagger: float = 0.0,
    highlight: Sequence | None = None,
    highlight_size_boost: float = 2.0,
    fallback_color: str = "lightgrey",
    add_outline: bool = False,
    outline_width: tuple[float, float] = (0.2, 0.025),
    outline_color: tuple[str, str] = ("black", "white"),
) -> str:
    """Render a video morphing cells from ``adata.obsm[obsm_from]`` to ``adata.obsm[obsm_to]``.

    The motivating case is spatial transcriptomics -> UMAP (``obsm_from="spatial"``,
    ``obsm_to="X_umap"``): watching each cell fly from where it physically sits
    in the tissue to where it lands transcriptomically. Any other ``.obsm`` pair
    (PCA -> UMAP, t-SNE -> UMAP, before/after batch correction, ...) works the
    same way -- just pass different ``obsm_from``/``obsm_to`` keys.

    Cells are matched by row order (both obsm arrays describe the same cells
    in the same AnnData, so no assignment/matching problem needs solving).
    When ``align`` is True, ``obsm_from`` is first Procrustes-aligned onto
    ``obsm_to``'s frame (rotation/reflection/scale) so cells flow coherently
    instead of scattering across the canvas -- recommended whenever the two
    embeddings aren't already in comparable orientation/scale (e.g. spatial
    vs UMAP, or PCA vs UMAP).

    Parameters
    ----------
    adata : AnnData
    obsm_from, obsm_to : keys into ``adata.obsm``
    color : key into ``adata.obs`` used to color cells (categorical or continuous)
    out_path : output video path, e.g. ``"transition.mp4"``
    n_frames : number of interpolated frames between the two embeddings
    fps : output frames per second
    hold_frames : frames to hold still on each end (lets a viewer read the plot).
        Defaults (``n_frames=130``, ``hold_frames=40``, ``fps=30``) add up to
        a 7-second video.
    align : Procrustes-align ``obsm_from`` onto ``obsm_to`` before interpolating
    dims : which two columns of each obsm array to use as (x, y)
    labels : display names for (obsm_from, obsm_to); defaults to the keys themselves
    palette : ``{category: color}`` mapping, a palette name (see
        ``embedmorph.palettes.list_palettes()``), a list of colors, or ``None``
        for the default. Ignored for continuous ``color`` columns.
    zoom_to_fit : frame the view tightly around whatever's actually on screen
        at each frame -- each hold sits tight around that embedding alone,
        and during the morph the camera continuously tracks the moving
        point cloud's own bounding box frame by frame (naturally widening
        only as far as cells are actually spread out at that instant, then
        settling back down as they converge on their destination). This is
        what makes each embedding fill the frame instead of sitting tiny in
        a canvas sized to fit both at once. Set False for one fixed view
        (spanning both embeddings) for the whole video instead.
    view_padding : fraction of each view's own coordinate range added as
        margin around it.
    arc_strength : cells travel along a gentle curve rather than a straight
        line, bulging perpendicular to their own direction of travel by this
        fraction of the distance covered (every cell curves the same way, so
        the cloud reads as one coherent flow). 0 (default) for straight-line
        travel; try ~0.15 for a noticeable "flock" curve.
    stagger : spreads cell launch times over this fraction of the morph's
        frames instead of every cell starting and stopping in perfect
        lockstep -- cells with farther to travel leave sooner, so everyone
        still lands together by the last frame. 0 (default) for a fully
        synchronized launch; try ~0.3 for an organic staggered departure.
    highlight : subset of ``color``'s categories to keep in color; every
        other cell is drawn ``fallback_color`` and dropped from the legend
        (categorical ``color`` only) -- e.g. to spotlight one cell type.
    highlight_size_boost : multiplier on ``point_size`` for highlighted cells.
    add_outline : draw each point with a thin ring (as scanpy's
        ``add_outline`` does), so cells stand out from the background and
        from each other -- costs roughly 3x the per-frame render time.
    outline_width, outline_color : as in scanpy's ``sc.pl.embedding``:
        ``outline_width`` is ``(ring_size, gap_size)`` as a fraction of the
        point radius, ``outline_color`` is ``(ring_color, gap_color)``.

    Returns
    -------
    str : the output path, for convenience.
    """
    if obsm_from not in adata.obsm:
        raise KeyError(f"{obsm_from!r} not found in adata.obsm; available: {list(adata.obsm.keys())}")
    if obsm_to not in adata.obsm:
        raise KeyError(f"{obsm_to!r} not found in adata.obsm; available: {list(adata.obsm.keys())}")
    if color not in adata.obs:
        raise KeyError(f"{color!r} not found in adata.obs; available: {list(adata.obs.columns)}")
    if easing not in _EASINGS:
        raise ValueError(f"unknown easing {easing!r}; choose from {list(_EASINGS)}")

    dim0, dim1 = dims
    coords_from = np.asarray(adata.obsm[obsm_from])[:, [dim0, dim1]].astype(float)
    coords_to = np.asarray(adata.obsm[obsm_to])[:, [dim0, dim1]].astype(float)
    if coords_from.shape[0] != adata.n_obs or coords_to.shape[0] != adata.n_obs:
        raise ValueError("obsm arrays must have one row per cell in adata")

    coords_from_aligned = procrustes_align(coords_from, coords_to) if align else coords_from
    color_mapping = resolve_colors(adata.obs[color], palette, highlight=highlight, fallback_color=fallback_color)

    view_from = _view_bounds(coords_from_aligned, view_padding)
    view_to = _view_bounds(coords_to, view_padding)

    label_from, label_to = labels or (obsm_from, obsm_to)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi, facecolor=facecolor)
    legend_visible = show_legend and color_mapping.kind == "categorical" and len(color_mapping.legend) <= 60
    right_margin = 0.72 if legend_visible else 0.97
    fig.subplots_adjust(left=0.03, right=right_margin, top=0.92, bottom=0.03)

    ease_fn = _EASINGS[easing]
    raw_us = np.linspace(0.0, 1.0, n_frames)

    control = _arc_controls(coords_from_aligned, coords_to, arc_strength) if arc_strength > 0 else None
    stagger_start = _stagger_starts(coords_from_aligned, coords_to, stagger)
    stagger_span = max(1.0 - stagger, 1e-9)

    def _positions_at(u: float) -> np.ndarray:
        local_u = np.clip((u - stagger_start) / stagger_span, 0.0, 1.0)
        t = ease_fn(local_u)[:, None]
        if control is None:
            return (1 - t) * coords_from_aligned + t * coords_to
        return (1 - t) ** 2 * coords_from_aligned + 2 * (1 - t) * t * control + t**2 * coords_to

    order, color_mapping, point_size = reorder_for_highlight(color_mapping, point_size, highlight_size_boost)

    def _draw(coords):
        return coords[order] if order is not None else coords

    out_path = str(out_path)
    writer = imageio.get_writer(out_path, fps=fps, codec="libx264", quality=8, macro_block_size=None)

    def frame_iter():
        if zoom_to_fit:
            for _ in range(hold_frames):
                yield _draw(coords_from_aligned), view_from, label_from
            for u in raw_us:
                pos = _positions_at(u)
                # Camera continuously tracks the moving cloud's own extent (tight
                # at t=0/1, naturally wider only while cells are actually spread
                # out mid-flight) instead of jumping to a fixed worst-case view.
                yield _draw(pos), _view_bounds(pos, view_padding), f"{label_from} → {label_to}"
            for _ in range(hold_frames):
                yield _draw(coords_to), view_to, label_to
        else:
            shared_view = _enclosing(
                _view_bounds(np.vstack([coords_from_aligned, coords_to]), view_padding), view_from, view_to
            )
            for _ in range(hold_frames):
                yield _draw(coords_from_aligned), shared_view, label_from
            for u in raw_us:
                yield _draw(_positions_at(u)), shared_view, f"{label_from} → {label_to}"
            for _ in range(hold_frames):
                yield _draw(coords_to), shared_view, label_to

    try:
        frames = list(frame_iter())
        for coords, (xlim, ylim), title in tqdm(frames, disable=not show_progress, desc=f"{obsm_from}->{obsm_to}"):
            rgb = render_frame(
                fig,
                ax,
                coords,
                color_mapping,
                xlim,
                ylim,
                point_size,
                alpha,
                title,
                facecolor,
                legend_visible,
                add_outline,
                outline_width,
                outline_color,
            )
            writer.append_data(rgb)
    finally:
        writer.close()
        plt.close(fig)

    return out_path
