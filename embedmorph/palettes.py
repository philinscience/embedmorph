"""Named color palettes for categorical coloring.

ColorBrewer (``Set1``, ``Set2``, ``Paired``, ``Dark2``), ``Tableau``, CARTOColors
``Bold``, and journal-style (``Cell``, ``Nature``, ``Science``, ``JAMA``) qualitative
palettes, ported verbatim from `scplotkit <https://github.com/philinscience/scplotkit>`_
so figures made with embedmorph and scplotkit can share a palette.
"""

from __future__ import annotations

import colorsys
from collections.abc import Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

NAMED_PALETTES: dict[str, list[str]] = {
    # -- ColorBrewer qualitative --
    "Set1": ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"],
    "Set2": ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"],
    "Paired": [
        "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
        "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
    ],
    "Dark2": ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"],
    # -- other qualitative --
    "Tableau": ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"],  # noqa: E501
    "Bold": ["#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695", "#CF1C90", "#F97B72"],  # noqa: E501
    # -- journal-inspired --
    "Cell": ["#C84C3A", "#2F7E8F", "#E1A22E", "#4E5A8A", "#5F9862", "#D07A6A", "#8B6FA8", "#7B8C9E", "#B85F7A", "#6B6B6B"],  # noqa: E501
    "Nature": ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"],  # noqa: E501
    "Science": ["#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"],  # noqa: E501
    "JAMA": ["#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"],
}


def list_palettes() -> list[str]:
    """Return the names of all built-in named palettes. See :func:`palette_colors`."""
    return list(NAMED_PALETTES)


def _extra_hues(n: int) -> list[tuple[float, float, float]]:
    """``n`` extra colors spread around the hue wheel via the golden-angle trick.

    Used to extend a named/explicit palette once its own colors run out,
    without pulling in seaborn just for its ``husl`` palette.
    """
    golden = 0.6180339887498949
    h = 0.0
    colors = []
    for _ in range(n):
        h = (h + golden) % 1.0
        colors.append(colorsys.hls_to_rgb(h, 0.55, 0.65))
    return colors


def _cmap_colors(name: str) -> list:
    cmap = plt.get_cmap(name)
    if hasattr(cmap, "colors"):  # ListedColormap, e.g. "tab20", "Set3"
        return list(cmap.colors)
    return [cmap(x) for x in [i / 19 for i in range(20)]]  # continuous cmap: sample 20 steps


def palette_colors(palette: str | Sequence | None, n: int) -> list:
    """Return at least ``n`` colors, as RGB(A) tuples, for a palette source.

    ``palette`` is one of: a name from :func:`list_palettes`, any matplotlib
    colormap name (``"tab20"``, ``"Set3"``, ``"viridis"``, ...), an explicit
    sequence of colors, or ``None`` for the default -- ``tab20`` chained with
    ``tab20b``/``tab20c`` (60 distinct colors before it needs to fall back
    further).
    """
    if palette is None:
        base = _cmap_colors("tab20") + _cmap_colors("tab20b") + _cmap_colors("tab20c")
    elif isinstance(palette, str):
        base = [to_rgb(c) for c in NAMED_PALETTES[palette]] if palette in NAMED_PALETTES else _cmap_colors(palette)
    else:
        base = [to_rgb(c) for c in palette]

    if len(base) >= n:
        return base
    return base + _extra_hues(n - len(base))
