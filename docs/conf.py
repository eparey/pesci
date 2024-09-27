# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Path setup --------------------------------------------------------------

# # If extensions (or modules to document with autodoc) are in another directory,
# # add these directories to sys.path here. If the directory is relative to the
# # documentation root, use os.path.abspath to make it absolute, like shown here.
# #
import os
import sys
sys.path.insert(0, os.path.abspath('../'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pesci'
copyright = '2024, Elise Parey'
author = 'Elise Parey'
release = 'pesci 0.1.0'

master_doc = 'index'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx_nefertiti", "sphinx-prompt", "sphinx.ext.autosectionlabel",
              "sphinx.ext.autodoc", 'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_mock_imports = ["matplotlib", 'seaborn', "sklearn", "numpy", "pandas",
                        "scipy", 'coloredlogs', 'datatable', 'tqdm', 'threadpoolctl',
                        'scanpy', 'networkx']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# import sphinx_nefertiti

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_nefertiti"

html_static_path = ['_static']

html_theme_options = {
    
    "doc_headers_font": "Exo",  # Default value, no need to provide.
}
