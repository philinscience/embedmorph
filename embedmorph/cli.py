"""Command-line entry point: embedmorph"""

from __future__ import annotations

import argparse

import anndata as ad

from .palettes import list_palettes
from .transition import make_transition_video


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="embedmorph",
        description="Render a video morphing cells from one AnnData .obsm embedding to another -- "
        "typically spatial coordinates into a UMAP embedding, though any .obsm pair is supported.",
    )
    p.add_argument("h5ad", help="Path to an .h5ad file")
    p.add_argument("--obsm-from", required=True, help="adata.obsm key to start from, e.g. spatial (or X_pca)")
    p.add_argument("--obsm-to", required=True, help="adata.obsm key to end at, e.g. X_umap")
    p.add_argument("--color", required=True, help="adata.obs column to color cells by, e.g. cell_type")
    p.add_argument("-o", "--out", default="transition.mp4", help="Output video path (default: transition.mp4)")
    p.add_argument(
        "--n-frames",
        type=int,
        default=130,
        help="Interpolated frames between the two embeddings (defaults, with --hold-frames and "
        "--fps, add up to a 7-second video)",
    )
    p.add_argument("--fps", type=int, default=30, help="Output frames per second")
    p.add_argument(
        "--hold-frames", type=int, default=40, help="Frames to hold still at each end"
    )
    p.add_argument("--no-align", action="store_true", help="Skip Procrustes alignment before interpolating")
    p.add_argument("--dims", type=int, nargs=2, default=(0, 1), metavar=("X", "Y"), help="obsm columns to use as x, y")
    p.add_argument("--point-size", type=float, default=8.0)
    p.add_argument("--alpha", type=float, default=0.85)
    p.add_argument("--figsize", type=float, nargs=2, default=(8.0, 7.0), metavar=("W", "H"))
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--facecolor", default="white")
    p.add_argument("--no-legend", action="store_true")
    p.add_argument("--easing", choices=["linear", "ease_in_out", "ease_smoother"], default="ease_smoother")
    p.add_argument("--label-from", default=None, help="Display name for --obsm-from (default: the key itself)")
    p.add_argument("--label-to", default=None, help="Display name for --obsm-to (default: the key itself)")
    p.add_argument("--quiet", action="store_true", help="Suppress the progress bar")
    p.add_argument(
        "--no-backed",
        action="store_true",
        help="Load the full .h5ad into memory instead of backed mode (backed mode reads "
        "obs/obsm eagerly but leaves .X on disk, which this tool never touches -- much "
        "faster/lighter for large files, on by default)",
    )
    p.add_argument(
        "--palette",
        default=None,
        help=f"Color palette for categorical `--color` columns: a built-in name "
        f"({', '.join(list_palettes())}), any matplotlib colormap name (tab20, Set3, ...), "
        f"or omit for the default.",
    )
    p.add_argument(
        "--no-zoom-to-fit",
        action="store_true",
        help="Keep one fixed camera view for the whole video instead of continuously tracking "
        "the moving point cloud's own extent (zoom-to-fit is what makes each embedding fill "
        "the frame instead of sitting tiny in a canvas sized for both; on by default)",
    )
    p.add_argument("--view-padding", type=float, default=0.08, help="Margin around each view, as a fraction of its range")
    p.add_argument(
        "--arc-strength",
        type=float,
        default=0.0,
        help="Cells travel along a gentle curve (bulging perpendicular to their own direction "
        "of travel by this fraction of the distance covered) instead of a straight line. 0 "
        "(default) for straight-line travel; try ~0.15 for a noticeable 'flock' curve.",
    )
    p.add_argument(
        "--stagger",
        type=float,
        default=0.0,
        help="Spread cell launch times over this fraction of the morph's frames instead of "
        "every cell starting/stopping in lockstep -- farther cells leave sooner, everyone "
        "lands together by the last frame. 0 (default) for a fully synchronized launch; try "
        "~0.3 for an organic staggered departure.",
    )
    p.add_argument(
        "--highlight",
        nargs="+",
        default=None,
        metavar="CATEGORY",
        help="Keep only these categories of --color in color; every other cell is greyed out "
        "and dropped from the legend (e.g. to spotlight one cell type)",
    )
    p.add_argument("--highlight-size-boost", type=float, default=2.0, help="--point-size multiplier for highlighted cells")
    p.add_argument("--fallback-color", default="lightgrey", help="Color for non-highlighted cells when --highlight is used")
    p.add_argument(
        "--add-outline",
        action="store_true",
        help="Draw a thin ring around each point (as scanpy's add_outline does) -- costs "
        "roughly 3x the per-frame render time",
    )
    p.add_argument("--outline-width", type=float, nargs=2, default=(0.2, 0.025), metavar=("RING", "GAP"))
    p.add_argument("--outline-color", nargs=2, default=("black", "white"), metavar=("RING", "GAP"))
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    adata = ad.read_h5ad(args.h5ad, backed=None if args.no_backed else "r")

    labels = None
    if args.label_from or args.label_to:
        labels = (args.label_from or args.obsm_from, args.label_to or args.obsm_to)

    out = make_transition_video(
        adata,
        obsm_from=args.obsm_from,
        obsm_to=args.obsm_to,
        color=args.color,
        out_path=args.out,
        n_frames=args.n_frames,
        fps=args.fps,
        hold_frames=args.hold_frames,
        align=not args.no_align,
        dims=tuple(args.dims),
        point_size=args.point_size,
        alpha=args.alpha,
        figsize=tuple(args.figsize),
        dpi=args.dpi,
        facecolor=args.facecolor,
        palette=args.palette,
        show_legend=not args.no_legend,
        easing=args.easing,
        labels=labels,
        show_progress=not args.quiet,
        zoom_to_fit=not args.no_zoom_to_fit,
        view_padding=args.view_padding,
        arc_strength=args.arc_strength,
        stagger=args.stagger,
        highlight=args.highlight,
        highlight_size_boost=args.highlight_size_boost,
        fallback_color=args.fallback_color,
        add_outline=args.add_outline,
        outline_width=tuple(args.outline_width),
        outline_color=tuple(args.outline_color),
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
