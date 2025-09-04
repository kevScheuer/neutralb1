# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "neutralb1"
copyright = "2025, Kevin Scheuer"
author = "Kevin Scheuer"

import neutralb1

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # for our Google style docstrings
    "sphinx.ext.viewcode",  # add links to highlighted source code
    "sphinx.ext.autosummary",  # automatically generate api
]

# ---Autosummary settings---
autosummary_generate = True
autosummary_imported_members = True

# ---Autodoc settings---
autodoc_typehints = "description"  # place typehints in the description of functions
# Don't show class signature with the class' name.
autodoc_class_signature = "separated"

# TODO: publish html to github pages

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
