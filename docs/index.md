# embedmorph

Morph a spatial transcriptomics scatter plot into its UMAP embedding,
colored by an `.obs` column (cell type, niche, cluster...), and render it as
an mp4. That's the motivating case: watching each cell fly from where it
physically sits in the tissue to where it lands transcriptomically is one of
the most intuitive ways to *show* what a cluster or niche call means, rather
than just asserting it in two side-by-side static plots. Built for
presentations: a "flying dots" transition between two ways of looking at the
same cells.

The same machinery works for any pair of AnnData `.obsm` embeddings, not
just spatial → UMAP -- PCA → UMAP, t-SNE → UMAP, before/after batch
correction, etc. are all fully supported, spatial → UMAP just happens to be
the case this package was built around.

Because both `.obsm` arrays describe the same cells in the same row order,
there's no cell-matching problem to solve -- cell `i` in one embedding *is*
cell `i` in the other. The real design problem is making the motion look
coherent rather than chaotic and framing each embedding so it actually fills
the video instead of sitting tiny in a canvas sized for both -- see
{doc}`installation` for how that's handled, and the demonstration notebooks
for what it looks like.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
demonstrations/01_quickstart
demonstrations/02_gallery
api
```

## At a glance

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

See the {doc}`demonstrations/01_quickstart` for a runnable end-to-end
example, the {doc}`demonstrations/02_gallery` for palettes, highlighting,
point outlines, and the optional flowing-arcs/staggered-launch flourish, or
the {doc}`api` for the complete function reference.
