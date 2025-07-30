.. _starting:

######################
Quickstart with XSNAP
######################

Prerequisites
===============

You’ll need to know a bit of `Python <https://www.python.org>`_, `NumPy <https://numpy.org>`_, 
`Matplotlib <https://matplotlib.org>`_, `Astropy <https://www.astropy.org>`_, 
and `XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_ (or `PyXspec <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/index.html>`_). 
Please refer to each of their own documentation for details.

Additionally, a little understanding on basic astronomical `photometry <https://en.wikipedia.org/wiki/Photometry_(astronomy)>`_, 
`spectroscopy <https://en.wikipedia.org/wiki/Astronomical_spectroscopy>`_, and `stellar physics <https://ads.harvard.edu/books/1989fsa..book/>`_ 
would help.

Installation
===============

Details are available on :doc:`install` but a quick command below can be used.

.. code-block:: bash

    pip install xsnap

.. note::
    
    It is best to use a new virtual environment, particularly `conda <https://anaconda.org/anaconda/conda>`_ environment, for XSNAP.

.. danger::

    You **must first** install `HEASoft + XSPEC <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/download.html>`_ and `PyXspec <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/buildinstall.html>`_

Importing XSNAP
=================

Fundamentally, you can import XSNAP by:

.. code-block:: python

    import xsnap

However, it can be tricky as it heavily relies on 
`PyXspec <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/index.html>`_, 
which is not available in `PyPI <https://pypi.org/>`_.

In this page, there will be two walkthroughs on importing XSNAP, i.e. through Jupyter and Visual Studio Code (VS Code)

Importing through Jupyter
---------------------------

1. Activate the virtual environment that you installed XSNAP in:

    >>> conda activate xsnap-venv

2. Initialize HEASoft in the environment:

    >>> export HEADAS=/path/to/heasoft/software
    >>> source $HEADAS/headas-init.sh

3. Run Jupyter lab or Jupyter Notebook

    >>> jupyter lab 

   or

    >>> jupyter notebook

4. Import XSNAP through the ``.ipynb`` file: 
    
    >>> import xsnap

Importing through VS Code
---------------------------

1. Make sure to add VS Code to the ``PATH`` environment variable, such that you can open VS Code using ``code`` in terminal. 
   Details can be found in the VS Code setup page for `Linux <https://code.visualstudio.com/docs/setup/linux#_install-vs-code-on-linux>`_, `macOS <https://code.visualstudio.com/docs/setup/mac#_launch-vs-code-from-the-command-line>`_, and `Windows <https://code.visualstudio.com/docs/setup/windows#_install-vs-code-on-windows>`_.

2. Activate the virtual environment that you installed XSNAP in:

    >>> conda activate xsnap-venv

3. Initialize HEASoft in the environment:

    >>> export HEADAS=/path/to/heasoft/software
    >>> source $HEADAS/headas-init.sh

4. Open VS Code through the exact same terminal:

    >>> code .

5. Import XSNAP through the ``.ipynb`` file: 
    
    >>> import xsnap