# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import astropy.units as u
from importlib.metadata import version

sys.path.insert(0, os.path.abspath("../src"))


project = 'XSNAP'
copyright = '2025, Ferdinand, Wynn V. Jacobson-Galan, and Mansi M. Kasliwal'
author = 'Ferdinand, Wynn V. Jacobson-Galan, Mansi M. Kasliwal'

# import xsnap

project = "XSNAP"
release = version("xsnap")
version = release



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
    'notfound.extension',
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
project = "XSNAP"
html_title = f"{project} v{version}"

html_static_path = ['_static']
html_css_files = ["custom.css"]
templates_path   = ["_templates"]

html_theme = "pydata_sphinx_theme"
html_show_sourcelink = False


html_theme_options = {
    "navbar_align": "content",
    "logo": {
        "image_light": "_static/xsnaplogo.svg",
        "image_dark": "_static/logo/xsnap_logo_header_transparent_white.png",
    },
    "github_url": "https://github.com/fercananything/xsnap",
    "external_links": [
        {"name": "Changelog", "url": "https://github.com/fercananything/XSNAP/blob/main/CHANGELOG.md"},
    ],
    "icon_links": [
        {
            "name": "E-mail",
            "url": "mailto:support@xsnap.org",
            "icon": "fas fa-inbox",   
            "type": "fontawesome",
        },
    ],
    "secondary_sidebar_items": ["page-toc"],
    "header_links_before_dropdown": 6,
    "collapse_navigation": True,
    "analytics": {
        "google_analytics_id": "G-MJMZJK59M0",  # GA4 measurement ID
    },
    # "show_toc_level": 2,
}

html_sidebars = {

}

notfound_context = {
    "title": "Page not found",
    "body": (
        "<h1>Oops! <br /> The page you are looking for is not found</h1>"
        "<h2>Please go back to the "
        "<a href='/index.html'>homepage</a>.</h2>"
    ),
}

# If you host at a subpath (e.g., GitHub Pages project site):
html_baseurl = "https://xsnap.org/"
notfound_urls_prefix = "/"

# ---------------- SEO extensions ----------------
extensions += [
    "sphinx_sitemap",          # sitemap.xml
    "sphinxext.opengraph",     # OG/Twitter cards
    "sphinxcontrib.robots",    # robots.txt
    "sphinxext.rediraffe",     # 301 redirects when pages move
]

# Required by sitemap & canonical logic (you already set this)
html_baseurl = "https://xsnap.org/"

# Sitemap: use pretty links that match how pages resolve
sitemap_url_scheme = "{link}"

# Robots: point to the sitemap; (optionally) disallow source dumps
robots_txt = f"""User-agent: *
Allow: /
Sitemap: {html_baseurl}sitemap.xml
Disallow: /_sources/
"""

# OpenGraph defaults (customize image if you have a social card)
ogp_site_url = html_baseurl
ogp_site_name = project
ogp_type = "website"
# Recommended: a 1200x630 PNG/JPG hosted on your site
ogp_image = "https://xsnap.org/_static/logo/xsnap_logo_icon_crop.jpeg"

# Treat broken xrefs as warnings (helps internal linking quality)
nitpicky = True

# If you ever rename/move pages, declare 301s here:
rediraffe_redirects = {
    # "old/path/index.rst": "new/path/",
    # "faq.rst": "guide/faq/",
}



