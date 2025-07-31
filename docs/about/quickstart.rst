.. _quickstart:

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

.. important::

    You **must remember to first** install `HEASoft + XSPEC <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/download.html>`_ and `PyXspec <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/buildinstall.html>`_.
    Details of the required dependencies are :ref:`here <required-dependencies>`.

    Before doing anything with XSNAP, **make sure** to initialize `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_:

    If you have Bourne shell variants (bash/sh/zsh):

        >>> export HEADAS=/path/to/heasoft/platform/
        >>> source $HEADAS/headas-init.sh

    Or if you have C-shell variants (tcsh/csh):

        >>> setenv HEADAS /path/to/heasoft/platform/
        >>> source $HEADAS/headas-init.csh

.. _hint-heasoft:
.. hint::
    
    It is best to use a new virtual environment, particularly `conda <https://anaconda.org/anaconda/conda>`_ environment, for XSNAP.

    In addition, it may be a good idea to make an alias for `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ 
    in your ``~/.bashrc`` or ``~/.zshrc`` by doing:

    #. Open the ``~/.bashrc`` or ``~/.zshrc`` file:
        
        >>> nano ~/.bashrc       # or ~/.zshrc

    #. Add these two lines:

        >>> export HEADAS=/path/to/heasoft/platform/
        >>> alias heainit='. $HEADAS/headas-init.sh'
    
    #. Reload the shell config:

        >>> source ~/.bashrc     # or ~/.zshrc

    #. Now, you can initialize `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ by:

        >>> heainit

    Please adjust accordingly if you have C-shell variants (tcsh/csh).
    

Quick CLI Usage
===========================

Once installed, XSNAP command-line interface (CLI) scripts are easy and ready-to-use. 

Initialize `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ in the terminal first:

    >>> export HEADAS=/path/to/heasoft/software
    >>> source $HEADAS/headas-init.sh

A quick :ref:`hint <hint-heasoft>` is also available above.

Then, you can immediately use in terminal, e.g.,

    >>> extract-xmm ./0882480901/odf/ sou_physical.reg bkg_physical.reg

If you want to run it in an ``.ipynb`` file:

    >>> ! extract-xmm ./0882480901/odf/ sou_physical.reg bkg_physical.reg

You can also run it in vanilla Python:

    >>> import subprocess
    >>> cmd = ["make-region", 
    ...        "/path/to/event/file", 
    ...        "169.59", "-32.83"]
    >>> subprocess.run(cmd)

.. note::

    Make sure to have install the :ref:`dependencies <recommended-dependencies>` of each script.
    Details of each script is also available :doc:`here <../api/index>`.
    Some examples are also made :doc:`here <../examples/notebook/Scripts>`.

Quick API Usage
========================

Importing XSNAP
^^^^^^^^^^^^^^^^^^^

Fundamentally, you can import XSNAP by:

.. code-block:: python

    import xsnap

However, it can be tricky as it heavily relies on 
`PyXspec <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/index.html>`_, 
which is not available in `PyPI <https://pypi.org/>`_.

In this page, there will be two walkthroughs on importing XSNAP, i.e. through `Jupyter <(tcsh/csh)>`_ and `Visual Studio Code (VS Code) <https://code.visualstudio.com>`_

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

Spectrum Fitting
^^^^^^^^^^^^^^^^^^^^^^^^

Below is a minimal example to fit a spectrum. 

.. hint::
    
    The :py:class:`~xsnap.SpectrumFit` assumes you have a grouped spectrum file ``myspectrum.pi`` or at least have the same name as the background and response file

.. code-block:: python

    import xsnap

    # 1) Create the fit object
    spec = xsnap.SpectrumFit(abund="aspl")

    # 2) Load your PHA spectrum file
    spec.load_data("myspectrum.pi")

    # 3) Define a model (e.g. absorbed power law)
    spec.set_model("tbabs*pow", 
                    TBabs_nH="0.05 -1", 
                    powerlaw_PhoIndex=2)
    # The model has 3 parameters here: 
    # TBabs.nH, powerlaw.PhoIndex, and powerlaw.norm

    # 4) Fit!
    spec.fit(nIterations=500)

    # 5) Plot (if you want to see the plot in svg for example)
    spec.set_plot("ldata", device="/svg")

    # 6) Getting best-fit parameters for powerlaw PhoIndex and norm.
    # This means that we're fitting for 1σ for parameter 2 and 3 (PhoIndex and norm)
    spec.get_params("1.0 2 3")

    # 7) Getting observation times and count rates
    spec.get_counts()
    spec.get_time()

    # 8) Getting flux and luminosity
    df_flux = spec.get_fluxes()
    spec.get_lumin(df_flux['unabsorbed'], redshift=0.1)

Handling a Collection of Spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is a minimal example to plot light curves and parameter evolutions from multiple spectra. 

.. code-block:: python

    import xsnap

    # 1) Create manager object
    manager = xsnap.SpectrumManager()

    # 2) Analyze a spectrum
    spec1 = xsnap.SpectrumFit(abund="aspl") 
    spec1.load_data("spectrum1.pha")
    spec1.set_model("tbabs*pow", 
                    TBabs_nH="0.05 -1", 
                    powerlaw_PhoIndex=2)
    spec1.fit()
    spec1.get_params("1.0 2 3")
    spec1.get_counts()
    spec1.get_time()
    df_flux = spec1.get_fluxes()
    spec1.get_lumin(df_flux['unabsorbed'], redshift=0.1)

    # 3) Load the spectrum to manager
    # (Optionally) you can also specify 
    # the instrument by manager.load(spec1, instrument="XMM")
    # Otherwise, the manager will try to find the instrument
    # from the spectrum file
    manager.load(spec1)  

    # 4) Repeat analysis of other spectrums
    # Due to the nature of PyXspec, it is
    # best to load spectrum each time after analysis

    spec2 = xsnap.SpectrumFit(abund="aspl") 
    # clear existing data before loading a new one by doing:
    spec2.load_data("spectrum2.pha", clear=True)
    # you can also clear by importing xspec and run:
    # xspec.AllData.clear()
    spec2.set_model("tbabs*pow", 
                    TBabs_nH="0.05 -1", 
                    powerlaw_PhoIndex=2)
    spec2.fit()
    spec2.get_params("1.0 2 3")
    spec2.get_counts()
    spec2.get_time()
    df_flux = spec2.get_fluxes()
    spec2.get_lumin(df_flux['unabsorbed'], redshift=0.1)

    # 5) Load the spectrum again to manager
    manager.load(spec2)

    # 6) Plot flux and luminosity light curves
    manager.plot_flux()
    manager.plot_lumin()

    # 7) Plot parameter evolution
    manager.plot_params()

Detecting Sources
^^^^^^^^^^^^^^^^^^^^^^^^

Below is a minimal example to detect a source given an event file and the source coordinates.

.. code-block:: python

    import glob
    from xsnap import SourceDetection

    # 1) Define the paths to the event files and exposure image files 
    # (exposure image are optional but helpful and recommended to use)
    # Here, I use glob to list all the .evt and .img file inside of a directory

    evt_paths = sorted(glob.glob("./parent/dir/of/event/*.evt"))
    img_paths = sorted(glob.glob("./parent/dir/of/image/*.img"))

    # Keep in mind that the amount of event path and image must be the same
    # i.e. len(evt_paths) == len(img_paths)
    # This is because they will be paired up when used for detecting, i.e.
    # evt_paths[0] with img_paths[0]

    # 2) Create the detection object and load the files
    detect = SourceDetection()
    detect.load(evt_paths, img_paths)

    # 3) Define the source RA and dec
    # Make sure it's in decimal degrees
    RA = 169.59
    dec = -32.83

    # 4) Detect the source in all event files!
    detect.detect_all(RA, dec)

More Examples
^^^^^^^^^^^^^^^^^^

More examples are available in the :doc:`../examples/index` page (more notebooks will be made soon!).