# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import astropy.units as u

sys.path.insert(0, os.path.abspath("../src"))


project = 'XSNAP'
copyright = '2025, Ferdinand, Wynn V. Jacobson-Galan, and Mansi M. Kasliwal'
author = 'Ferdinand, Wynn V. Jacobson-Galan, Mansi M. Kasliwal'

import xsnap

# Sphinx’s “short” version (major.minor)
version = xsnap.__version__
# Sphinx’s “full” release string (can include patch, alpha/beta tags)
release = xsnap.__version__

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_design",
    "sphinx.ext.intersphinx",
    "nbsphinx",
]

autodoc_mock_imports = [
    "numpy",
    "pandas",
    "astropy",
    "astropy.time",
    "astropy.units",
    "erfa",
    "matplotlib",
    "matplotlib.pyplot",
    "emcee",
    "corner",
    "scipy",
    "xspec"
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "scipy":  ("https://docs.scipy.org/doc/scipy/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
}

# autodoc_typehints = "none"

autodoc_default_options = {
  "members": True,
  "undoc-members": True,
  "show-inheritance": True,
}

add_module_names = False

napoleon_numpy_docstring = True
napoleon_google_docstring = False

napoleon_use_param = True
napoleon_use_rtype = False  

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_favicon = '_static/xsnapicon.ico'

html_static_path = ['_static']
html_css_files = ["custom.css"]
templates_path   = ["_templates"]

html_theme = "pydata_sphinx_theme"
html_show_sourcelink = False


html_theme_options = {
    "navbar_align": "content",
    "logo": {
        "image_light": "_static/xsnaplogo.svg",
        "image_dark": "_static/xsnaplogo_dark.svg",
    },
    "github_url": "https://github.com/fercananything/xsnap",
    "external_links": [
        {"name": "Changelog", "url": "https://github.com/fercananything/XSNAP/blob/main/CHANGELOG.md"},
    ],
    "secondary_sidebar_items": ["page-toc"],
    "header_links_before_dropdown": 6,
    "collapse_navigation": True,
    # "show_toc_level": 2,
}

html_sidebars = {

}



