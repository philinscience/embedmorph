# embedmorph

[![Docs](https://github.com/philinscience/embedmorph/actions/workflows/docs.yml/badge.svg)](https://embedmorph.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](pyproject.toml)

![embedmorph demo: cells morphing from their spatial coordinates in tissue to a UMAP embedding](docs/_static/embedmorph_demo.gif)

Morph a spatial transcriptomics scatter plot into its UMAP embedding, colored
by an `.obs` column (cell type, niche, cluster...), and render it as an mp4.
That's the motivating case: watching each cell fly from where it physically
sits in the tissue to where it lands transcriptomically is one of the most
intuitive ways to *show* what a cluster or niche call means, rather than
just asserting it in two side-by-side static plots. Built for presentations --
a "flying dots" transition between two ways of looking at the same cells.

The same machinery works for any pair of AnnData `.obsm` embeddings, not
just spatial → UMAP -- PCA → UMAP, t-SNE → UMAP, before/after batch
correction, etc. are all fully supported (just pass different `obsm_from`/
`obsm_to` keys), spatial → UMAP just happens to be the case this package was
built around.

Because both `.obsm` arrays describe the same cells in the same row order,
there's no cell-matching problem to solve -- cell `i` in one embedding *is*
cell `i` in the other. The only real design problem is making the motion
look coherent rather than chaotic, which is handled with:

- **Procrustes alignment**: before interpolating, the starting embedding is
  rotated/reflected/scaled onto the target embedding's frame (translation +
  rotation + isotropic scale, least-squares optimal). Without this, two
  embeddings that happen to be flipped or rotated relative to each other
  produce a video where the whole cloud swings across the canvas instead of
  cells flowing directly to their destination.
- **Eased interpolation**: a smootherstep curve (zero velocity *and* zero
  acceleration at both ends) instead of linear interpolation, so the
  transition doesn't start/stop abruptly.
- **Zoom-to-fit camera** (`zoom_to_fit=True`, default): each hold sits tight
  around that embedding alone, and during the morph the camera continuously
  tracks the *moving point cloud's own bounding box*, frame by frame --
  naturally widening only as far as cells are actually spread out at that
  instant, then settling back down as they converge on their destination,
  rather than jumping to one fixed worst-case view for the whole morph.
  Without this, an embedding whose Procrustes-aligned scale happens to be
  small relative to the other (common for spatial vs UMAP -- a cell's
  physical location and its transcriptomic UMAP position are only weakly
  correlated, so the least-squares alignment scale ends up small) renders as
  a tiny cluster in a canvas sized for the bigger one. Set `zoom_to_fit=False`
  to keep one fixed view (spanning both embeddings) for the whole video instead.

Also built in, all off by default so the base look stays clean and
predictable: **flowing arcs + staggered launch** (`arc_strength=0.15`,
`stagger=0.3` are good starting points) for a more organic "flock" motion --
cells travel along a gentle curve rather than a straight line (every cell
bulges to the same side of its own direction of travel, so the cloud reads
as one coherent flow), and don't all launch in lockstep, though everyone
still lands together by the last frame; **more color palettes** than
matplotlib's defaults -- ColorBrewer, Tableau, and journal-style (`Nature`,
`Science`, `Cell`, `JAMA`) qualitative palettes, ported from
[scplotkit](https://github.com/philinscience/scplotkit)
(`embedmorph.list_palettes()` for the full list, pass a name as `palette=`);
**highlighting** one or more categories against a greyed-out background
(`highlight=[...]`, as in scplotkit's `masked_umap_highlight`); and
scanpy-style **point outlines** (`add_outline=True`, same construction as
`sc.pl.embedding`'s own `add_outline`) so cells stand out from a dense
background.

## Install

```bash
pip install -e .
# or, for the tqdm progress bar during rendering:
pip install -e ".[progress]"
```

Requires `ffmpeg` to be available to `imageio-ffmpeg` (it bundles its own
static binary by default, so this normally works out of the box even on a
compute node with no system ffmpeg module).

## Usage

### As a library

```python
from embedmorph import make_transition_video

make_transition_video(
    adata,
    obsm_from="spatial",
    obsm_to="X_umap",
    color="cell_type",
    out_path="transition.mp4",
)
```

(Any other `.obsm` pair -- `obsm_from="X_pca"`, etc. -- works exactly the
same way.) See `docs/demonstrations/01_quickstart.ipynb` for a runnable
end-to-end example (synthetic spatial data, plus a snippet for real data via
`scanpy`), and `docs/demonstrations/02_gallery.ipynb` for palettes,
highlighting, point outlines, and the flowing-arcs/staggered-launch
flourish.

### From the command line

```bash
embedmorph my_spatial_data.h5ad \
  --obsm-from spatial --obsm-to X_umap \
  --color cell_type \
  -o transition.mp4
```

Run `embedmorph --help` for all options (frame count, fps, point size,
figure size/dpi, easing, legend, colors, hold duration, etc).

A few worth calling out:

```bash
# Journal-style palette
embedmorph my_data.h5ad --obsm-from spatial --obsm-to X_umap \
  --color cell_type -o transition.mp4 --palette Nature

# Spotlight one or two cell types against a greyed-out atlas
embedmorph my_data.h5ad --obsm-from spatial --obsm-to X_umap \
  --color cell_type -o transition.mp4 \
  --highlight "CD8 T cells" "Tregs" --highlight-size-boost 2.5

# Organic flock motion instead of plain straight-line/synchronized travel
embedmorph my_data.h5ad --obsm-from spatial --obsm-to X_umap \
  --color cell_type -o transition.mp4 \
  --arc-strength 0.15 --stagger 0.3

# scanpy-style point outlines (costs ~3x the per-frame render time --
# best for smaller cell counts or a subsample, gets muddy on dense atlases)
embedmorph my_data.h5ad --obsm-from spatial --obsm-to X_umap \
  --color cell_type -o transition.mp4 --add-outline
```

### On a compute cluster

Rendering is single-process (no GPU needed); the main memory driver is the
number of cells the backed AnnData holds in `.obs`/`.obsm`, not `.X` --
`embedmorph` never touches `.X`, so the CLI reads the file in backed mode by
default. For anything past a few thousand cells, render on a compute node
rather than interactively; if you're rendering many datasets, submit one job
per video rather than looping inside one job.

## Documentation

Full API reference and demonstration notebooks live in
[`docs/`](docs/) and are published at https://embedmorph.readthedocs.io. To
build the docs locally:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs docs/_build/html
```

## Notes / tradeoffs

- Currently uses the first two dimensions of each `.obsm` array (configurable
  via `--dims`); this is a 2D scatter transition, not a 3D one.
- `align=True` (default) is recommended whenever the two embeddings aren't
  already in a comparable orientation/scale -- which is the common case
  whenever a UMAP is on one side, whether the other side is spatial
  coordinates or PCA/t-SNE. Turn it off (`--no-align`) if you specifically
  want to see the raw, unaligned coordinate difference.
- Legends are only drawn for categorical colors with <=60 categories (more
  than that stops being readable in a corner legend). The default palette
  chains `tab20`+`tab20b`+`tab20c` (60 distinct colors) before it starts
  repeating.

## License

MIT -- see [LICENSE](LICENSE).
