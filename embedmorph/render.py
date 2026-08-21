"""Color resolution and single-frame rendering for transition videos."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence

import matplotlib

matplotlib.use("Agg")  # headless: must happen before pyplot is imported

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

from .palettes import palette_colors


@dataclass
class ColorMapping:
    kind: str  # "categorical" | "continuous"
    colors: np.ndarray  # (n_cells, 4) RGBA
    legend: object  # list[(label, rgba)] for categorical, (cmap, norm) for continuous
    highlight_mask: np.ndarray | None = None  # (n_cells,) bool, only set when `highlight` was used


def resolve_colors(
    values: pd.Series,
    palette: dict | str | Sequence[str] | None = None,
    highlight: Sequence | None = None,
    fallback_color: str = "lightgrey",
) -> ColorMapping:
    """Map a ``.obs`` column to per-cell RGBA colors, categorical or continuous.

    ``palette`` is a ``{category: color}`` mapping (colors pinned per
    category, missing ones auto-filled), the name of a built-in or matplotlib
    palette (see :func:`embedmorph.palettes.list_palettes`), an explicit list of
    colors, or ``None`` for the default. Only meaningful for categorical
    columns -- continuous columns always use ``viridis``.

    ``highlight``, if given, keeps only those category values in color;
    every other cell is drawn ``fallback_color`` and dropped from the legend
    -- e.g. to spotlight one cell type against the rest of the atlas.
    """
    is_numeric = pd.api.types.is_numeric_dtype(values) and not (
        pd.api.types.is_bool_dtype(values) or isinstance(values.dtype, pd.CategoricalDtype)
    )

    if is_numeric:
        if highlight is not None:
            raise ValueError("`highlight` only applies to categorical `.obs` columns, not continuous ones.")
        arr = values.to_numpy(dtype=float)
        vmin, vmax = np.nanmin(arr), np.nanmax(arr)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("viridis")
        colors = cmap(norm(arr))
        return ColorMapping(kind="continuous", colors=colors, legend=(cmap, norm))

    categories = pd.Categorical(values)
    cats = list(categories.categories)

    if isinstance(palette, dict):
        missing = [c for c in cats if c not in palette]
        fill = dict(zip(missing, palette_colors(None, len(missing)))) if missing else {}
        cat_colors = [to_rgba(palette[c]) if c in palette else to_rgba(fill[c]) for c in cats]
    else:
        base = palette_colors(palette, len(cats))
        cat_colors = [to_rgba(base[i]) for i in range(len(cats))]

    highlight_mask = None
    legend_cats = cats
    if highlight is not None:
        highlight_set = set(highlight)
        unknown = [h for h in highlight_set if h not in cats]
        if unknown:
            raise ValueError(f"`highlight` values not found in the column's categories: {unknown}")
        fallback_rgba = to_rgba(fallback_color)
        cat_colors = [c if cat in highlight_set else fallback_rgba for cat, c in zip(cats, cat_colors)]
        legend_cats = [c for c in cats if c in highlight_set]
        highlight_mask = values.isin(highlight_set).to_numpy()

    unknown_color = to_rgba((0.7, 0.7, 0.7))
    colors = np.array([cat_colors[c] if c >= 0 else unknown_color for c in categories.codes])
    legend = [(str(c), cat_colors[cats.index(c)]) for c in legend_cats]
    return ColorMapping(kind="categorical", colors=colors, legend=legend, highlight_mask=highlight_mask)


def render_frame(
    fig: "plt.Figure",
    ax: "plt.Axes",
    coords: np.ndarray,
    color_mapping: ColorMapping,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    point_size,
    alpha: float,
    title: str,
    facecolor: str,
    show_legend: bool,
    add_outline: bool = False,
    outline_width: tuple[float, float] = (0.2, 0.025),
    outline_color: tuple[str, str] = ("black", "white"),
) -> np.ndarray:
    """Draw one frame and return it as an (H, W, 3) uint8 RGB array.

    ``point_size`` may be a scalar or a per-cell array (already in ``coords``'
    row order). ``add_outline`` draws each point with a thin ring, the same
    construction as scanpy's own ``add_outline``: two flat-colored rings
    grown outward from the point radius by ``outline_width`` (fractions of
    it), drawn behind the point -- makes cells stand out from the background
    and from each other, especially on dense/overlapping regions.
    """
    ax.clear()
    ax.set_facecolor(facecolor)

    if add_outline:
        bg_width, gap_width = outline_width
        bg_color, gap_color = outline_color
        point_radius = np.sqrt(np.asarray(point_size, dtype=float))
        bg_size = (point_radius + point_radius * bg_width * 2) ** 2
        gap_size = (point_radius + point_radius * gap_width * 2) ** 2
        ax.scatter(coords[:, 0], coords[:, 1], s=bg_size, c=bg_color, linewidths=0)
        ax.scatter(coords[:, 0], coords[:, 1], s=gap_size, c=gap_color, linewidths=0)

    ax.scatter(coords[:, 0], coords[:, 1], s=point_size, c=color_mapping.colors, alpha=alpha, linewidths=0)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    text_color = "0.9" if facecolor in ("black", "#000000") else "0.1"
    if title:
        ax.set_title(title, color=text_color, fontsize=13, fontweight="bold")

    if show_legend and color_mapping.kind == "categorical":
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="", color=rgba, label=label, markersize=6)
            for label, rgba in color_mapping.legend
        ]
        ax.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=8,
            labelcolor=text_color,
        )

    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return buf[:, :, :3].copy()


def reorder_for_highlight(color_mapping: ColorMapping, point_size, highlight_size_boost: float):
    """Sort cells so highlighted ones draw last (on top of the greyed-out background).

    Returns ``(order, color_mapping, point_size)`` with ``color_mapping.colors``
    and ``point_size`` already permuted by ``order`` -- apply the same
    ``order`` to each frame's coordinates before calling :func:`render_frame`.
    If ``color_mapping.highlight_mask`` is ``None`` (no highlighting in use),
    returns everything unchanged and ``order=None``.
    """
    mask = color_mapping.highlight_mask
    if mask is None:
        return None, color_mapping, point_size

    order = np.argsort(mask.astype(int), kind="stable")  # non-highlighted first, highlighted (drawn on top) last
    sizes = np.where(mask, np.asarray(point_size, dtype=float) * highlight_size_boost, point_size)
    return order, replace(color_mapping, colors=color_mapping.colors[order]), sizes[order]
