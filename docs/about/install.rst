.. _install:


###############
Installation
###############

Installing XSNAP
==================

We strongly recommend that you make use of Python virtual environments, 
or (even better) Conda virtual environments when installing XSNAP.

Currently, `XSNAP <https://pypi.org/p/xsnap>`_ is available for download in the popular Python Package Index (PyPI) .

.. code-block:: bash

   pip install xsnap

Additionally, XSNAP should be able to be downloaded by cloning this Github repository and run:

.. code-block:: bash

   git clone https://github.com/fercananything/XSNAP/
   cd XNAP
   python -m pip install .

.. _dependencies:
Dependencies
===============

.. _required-dependencies:
Required Dependencies
------------------------

XSNAP analysis depends heavily on two non-Python softwares:

* `HEASOFT <https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html>`_ - Version 6.35. 
   Other recent versions should be compatible even if I have yet to test it. It is best to install all the packages; otherwise, install at least ``NuSTAR``, ``Swift``, all of ``General-Use FTOOLS``, and all of ``XANADU``.
* `HEASOFT's PyXspec <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/buildinstall.html>`_ - Version 2.1.4 (or XSPEC - Version 12.15.0). 
   Other recent versions should be compatible even if I have yet to test it. 
   Additionally, PyXspec should be automatically installed when you install HEASOFT.

.. _recommended-dependencies:
Recommended Dependencies
-------------------------

While it's not necessarily required, it is strongly recommended to download these non-Python softwares:

* `Chandra Interactive Analysis of Observations (CIAO) <https://cxc.harvard.edu/ciao/download/index.html>`_ - Version 4.17. 
   CIAO is needed if you want to do the spectral extraction from CXO data (running ``extract-chandra``). 
   It is recommended to install CIAO using the ``conda create`` command, i.e. install on a different 
   Python/Conda virtual environment. This is to seperate HEASOFT (and XSPEC) with CIAO and avoid clashes between modules. 
* `XMM Science Analysis System (SAS) <https://www.cosmos.esa.int/web/xmm-newton/sas-download>`_ - Version 22.1. 
   However, other recent versions should still be compatible. SAS is needed if you want to do data calibration and 
   spectral extraction for XMM-Newton. A few extra steps for SAS installation can be found 
   `here <https://www.cosmos.esa.int/web/xmm-newton/sas-thread-startup#>`_. This is needed to run ``extract-xmm``.
* `HEASARC Calibration Database (CALDB) <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/install.html>`_ - Version 2009 Aug 04. 
   The HEASARC CALDB is needed if you want to do data calibration and spectral extraction for Swift-XRT and NuSTAR.
   This is needed to run ``extract-nustar``, ``extract-swift``, and ``swift-stack-pc``.
* `CALDB Files for Swift-XRT and NuSTAR <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_supported_missions.html>`_. 
   In addition to the CALDB, the CALDB files are needed to be downloaded too. These files are needed if you want to do 
   data calibration and spectral extraction for Swift-XRT and NuSTAR. This is needed to run ``extract-nustar`` and ``extract-swift``.

.. caution::
   Keep in mind, without these softwares, you are only able to import the spectra fitting and analysis modules. 
   These softwares help with the scripts dealing for data calibration and spectral extraction.

.. _optional-dependencies:
Optional Dependencies
------------------------

This software is completely optional and has minimal impact on the user experience.

* `DS9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`_ - Version 4.1 and above. 
   DS9 is needed to help user's interactivity in making region files. This is (optionally) useful to run ``make-region``.

