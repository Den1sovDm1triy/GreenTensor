# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Импортируемый пакет лежит в корне репозитория (на уровень выше docs/).
# The import package lives at the repository root (one level above docs/).
sys.path.insert(0, os.path.abspath(".."))

import green_tensor  # noqa: E402

# -- Project information ------------------------------------------------------

project = "GreenTensor"
copyright = "2025-2026, GreenTensor authors (Ural Federal University)"
author = "D.V. Denisov, V.Ya. Noskov, I.O. Skumatenko"

version = green_tensor.__version__
release = green_tensor.__version__

# -- General configuration ----------------------------------------------------

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output --------------------------------------------------

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output --------------------------------------------------

epub_show_urls = "footnote"
