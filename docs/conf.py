from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import embedmorph  # noqa: E402

project = "embedmorph"
copyright = "2026, embedmorph contributors"
author = "embedmorph contributors"
release = embedmorph.__version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",
    "myst_nb",
]

templates_path = ["_templates"]
html_static_path = ["_static"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": False,
}
autodoc_typehints = "description"
napoleon_google_docstring = False
napoleon_numpy_docstring = True

nb_execution_mode = "off"  # demonstration notebooks are committed pre-executed
myst_enable_extensions = ["colon_fence"]

html_theme = "sphinx_book_theme"
html_title = "embedmorph"
html_css_files = ["custom.css"]
html_theme_options = {
    "repository_url": "https://github.com/philinscience/embedmorph",
    "use_repository_button": True,
    "use_download_button": False,
}
