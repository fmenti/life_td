# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
import os

# import pathlib
import sys

# sys.path.insert(0, pathlib.Path(__file__).parents[1].resolve().as_posix())
sys.path.insert(0, os.path.abspath("../../life_td_data_generation"))
sys.path.insert(0, os.path.abspath("tutorials"))


# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "life_td_data_generation"
copyright = "2024, Franziska Menti"
author = "Franziska Menti"

# The short X.Y version
with open(
    os.path.join(
        os.path.dirname(__file__), "../../life_td_data_generation/VERSION"
    )
) as version_file:
    version = version_file.read().strip()

# The full version, including alpha/beta/rc tags
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # Core library for html generation from docstrings
    "sphinx.ext.autosummary",  # Create neat summary tables
    "nbsphinx",
]  # jupyter notebooks

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"  #'alabaster'
html_static_path = ["_static"]
