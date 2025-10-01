# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the project root directory to Python path so scripts can be imported
sys.path.insert(0, os.path.abspath(".."))

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
    "sphinx.ext.mathjax",  # render math
    "sphinx.ext.todo",  # support for todo directives
]

# ---Autosummary settings---
autosummary_generate = True
autosummary_imported_members = False  # Don't include imported members

# ---Autodoc settings---
autodoc_typehints = "description"  # place typehints in the description of functions
# Don't show class signature with the class' name.
autodoc_class_signature = "separated"
# Include all members by default, including methods
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
# Skip imported members from other modules
autodoc_mock_imports = []

# Suppress duplicate warnings by not including members in module pages when they have their own pages
suppress_warnings = ["autosummary.import_cycle"]

# ---Todo extension settings---
todo_include_todos = True  # Include todo items in the build

# TODO: publish html to github pages

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
