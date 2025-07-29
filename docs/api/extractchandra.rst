.. _extract-chandra:

**********************************************
Extract CXO Data (:code:`extract-chandra`)
**********************************************

Pipeline module to generate :math:`500–8000` eV images, PSF maps, and 
extracted spectra from a Chandra X-ray Observatory (CXO) event file.

Automates the standard `CIAO <https://cxc.harvard.edu/ciao/>`_ workflow:
  1. Locate the user’s Conda environment.
  2. Locate the appropriate aspect solution (``.asol1.fits``) and mask (``.msk1.fits``) files.
  3. Run ``dmcopy`` to filter CCD 7 events in :math:`500–8000` eV and bin to a smaller image.
  4. Create a PSF map at 1.4967 keV with 90% encircled-energy fraction.
  5. Run ``specextract`` to build source and background spectra (NUM_CTS grouping), applying correct weights.

Positional arguments:
  - ``COND_ENV`` :  Name of the Conda environment to check against ``$CONDA_DEFAULT_ENV``.
  - ``EVT_FILE`` : Path to the Chandra event file.
  - ``SRC_REG`` :   Source region (filename or literal sky region).
  - ``BKG_REG`` :   Background region (filename or literal sky region).
  - ``OUTDIR`` :    Optional directory for outputs (defaults to ``EVT_FILE`` directory).

Options:
  --help          Show help message

Output files (inside ``OUTDIR``):
  - ``spec{obsid}_bkg.arf``
  - ``spec{obsid}_bkg.rmf``
  - ``spec{obsid}_bkg.pi``
  - ``spec{obsid}.arf``
  - ``spec{obsid}.rmf``
  - ``spec{obsid}.pi``
  - ``spec{obsid}.corr.arf``
  - ``spec{obsid}_grp.pi`` - the spectrum file that will be used in :py:class:`~xsnap.spectrum.SpectrumFit`

.. note::

    Assuming that the event file is inside ``obsid/primary/event.fits``, 
    the obsid is fundamentally the parent of the parent directory of the event file.

Usage formula:
    .. code-block:: bash

        extract-chandra <COND_ENV> <EVT_FILE> <SRC_REG> <BKG_REG> [OUTDIR]

Usage examples:
    .. code-block:: bash
        
        extract-chandra myenv /data/obs/1234/primary/event.fits src.reg bkg.reg
        extract-chandra myenv /data/obs/1234/primary/event.fits src.reg bkg.reg ./output 
        extract-chandra -h # for help

Requirements:
  • `CIAO <https://cxc.harvard.edu/ciao/>`_ (including ``dmcopy``, ``mkpsfmap``, ``specextract``) installed.
  • ``CONDA_DEFAULT_ENV`` set to the requested environment name.
