# Installation

## From source

```bash
git clone https://github.com/philinscience/embedmorph.git
cd embedmorph
pip install -e .

# or, for the tqdm progress bar during rendering:
pip install -e ".[progress]"
```

Requires `ffmpeg` to be available to `imageio-ffmpeg` (it bundles its own
static binary by default, so this normally works out of the box even on a
compute node with no system ffmpeg module).

## Requirements

embedmorph requires Python 3.9+ and works with any recent `anndata`
installation. It doesn't require GPU support or any compiled dependencies.

## Verifying the install

```python
import embedmorph
print(embedmorph.__version__)
```

## Command line

Installing the package also installs an `embedmorph` command:

```bash
embedmorph my_spatial_data.h5ad \
  --obsm-from spatial --obsm-to X_umap \
  --color cell_type \
  -o transition.mp4
```

Run `embedmorph --help` for the full option list (frame count, fps, point
size, figure size/dpi, easing, legend, colors, hold duration, palettes,
highlighting, point outlines, flowing arcs, and more).

## On a compute cluster

Anything past a few thousand cells is worth rendering on a compute node
rather than interactively. Rendering is single-process (no GPU needed); the
main memory driver is the number of cells the backed AnnData holds in
`.obs`/`.obsm`, not `.X` -- `embedmorph` never touches `.X`, so the CLI
reads the file in backed mode (`--no-backed` to opt out) rather than loading
the full expression matrix. If you're rendering many datasets, submit one
job per video rather than looping inside one job.
