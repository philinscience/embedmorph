# embedmorph

Morph a spatial transcriptomics scatter plot into its UMAP embedding,
colored by an `.obs` column (cell type, niche, cluster...), and render it as
an mp4. 

The same machinery works for any pair of AnnData `.obsm` embeddings, not
just spatial → UMAP -- PCA → UMAP, t-SNE → UMAP, before/after batch
correction, etc. are all fully supported, spatial → UMAP just happens to be
the case this package was built around.

see {doc}`installation` on how to install

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
