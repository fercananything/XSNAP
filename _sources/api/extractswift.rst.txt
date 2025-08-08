.. _extract-swift::

**********************************************
Extract Swift/XRT Data (:code:`extract-swift`)
**********************************************

Pipeline module to extract :math:`0.3–10.0` keV Swift/XRT spectra in PC (Photon Counting) or WT (Window Timing) mode.

Automates the standard `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ workflow for Swift/XRT:
  1. (Optionally) run `xrtpipeline <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/help/xrtpipeline.html>`_ to calibrate a raw data. (It is optional because it may have been run in `swift-stack-pc <swiftstackpc.html>`_)
  2. Filter cleaned events with user-supplied source / background regions.
  3. Create exposure-corrected ARFs and grouped PHA files.

Positional arguments:
  - ``OBSID``   :    11-digit Swift observation ID.
  - ``SRC_REG``  :    Source region (DS9 “physical” or “sky”).
  - ``BKG_REG``   :   Background region.

  .. hint::

    User can get source and background regions by running ``make-region`` in the clean event file (``*_cl.evt``) inside
    ``obsid/xrt/event/`` directory downloaded from raw Swift/XRT data. 
    Thus, there is no need to run ``xrtpipeline`` to get the region files.

    By default and recommended that user runs ``extract-swift`` on the parent of ``"sources"`` directory 
    and the ``"sources"`` directory has the downloaded data with ``ObsId`` as its child directory. E.g.

    .. code-block:: bash

        parent
        ├── sources
        │   ├── ObsId1
        │   └── ObsId2
        

    When ``extract-swift`` is run inside parent directory, it will result in

    .. code-block:: bash

        parent
        ├── sources
        │   ├── ObsId1
        │   └── ObsId2
        ├── output
        │   ├── ObsId1
        │   └── ObsId2
        

Options:
  --indir            Directory with event and exposure image files (default: ``./<OBSID>``).
  --outdir           Destination for spectrum file products (default: ``./`` if ``--no-pipe`` or ``./output/<OBSID>``).
  --src-ra           Source right ascension in ICRS (if not parsed, will try to use the RA in the source region file)
  --src-dec          Source declination in ICRS (if not parsed, will try to use the RA in the source region file)
  --mode             Extraction mode (PC default, options: PC or WT).
  --no-pipe          Bypass ``xrtpipeline``; use existing files.
  --evt              Path to event file
  --img              Path to exposure image file
  --help             Show help message

Usage formula:
    >>> extract-swift OBSID SRC_REG BKG_REG [--indir DIR] [--outdir DIR] \
    ...     [--src-ra RA] [--src-dec DEC] [--mode {PC | WT}] [--no-pipe] \
    ...     [--evt /path/to/evt] [--img /path/to/img]

Usage examples:
    Full calibration + extraction (PC mode):

      >>> extract 00012345001 src.reg bkg.reg --src-ra 150.123 --src-dec -12.345 --mode PC

    Skip the pipeline, work on existing events/images:

      >>> extract --no-pipe 00012345001 src.reg bkg.reg \
      ...     --indir ./events --evt sw00012345001xpcw3po_cl.evt

Output files:
    1. PC mode: 
        - ``pcsou.pha``  
        - ``pcbkg.pha``  
        - ``pcsou.arf``  
        - ``pcsougr1.pha`` - the spectrum file that will be used in :py:class:`~xsnap.spectrum.SpectrumFit`
    2. WT mode: 
        - ``wtsou.pha``  
        - ``wtbkg.pha``  
        - ``wtsou.arf``  
        - ``wtsougr1.pha`` - the spectrum file that will be used in :py:class:`~xsnap.spectrum.SpectrumFit`

Requirements:
  - `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ with `Swift/XRT CALDB files <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_supported_missions.html>`_ installed