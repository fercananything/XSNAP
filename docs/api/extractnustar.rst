.. _extract-nustar::

**********************************************
Extract NuSTAR Data (:code:`extract-nustar`)
**********************************************

Pipeline module to run NuSTAR reduction and spectral extraction in one step.

Automates the standard `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ workflow for NuStar:
  1. Execute `nupipeline <https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/nupipeline.html>`_ for a given ``OBSID``, producing cleaned ``FPMA`` and ``FPMB`` events.
  2. Run `nuproducts <https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/nuproducts.html>`_ to generate source and background spectra plus response files (ARFs/RMFs).
  3. Optionally reuse pre-cleaned events (skip pipeline).

Positional arguments:
  - ``OBSID``  :   11-digit NuSTAR observation ID.
  - ``SRC_REG`` :  Source region file or literal DS9 region string.
  - ``BKG_REG`` :  Background region file or literal DS9 region string.

Options:
  --indir         Directory containing raw event files (default: ``./sources/<OBSID>``).
  --outdir        Output root for products (default: ``./products/<OBSID>``).
  --ra            Source right ascension (RA) coordinates in ICRS (decimal degrees)
  --dec           Source declination (dec) coordinates in ICRS (decimal degrees)
  --no-pipe       Skip nupipeline; assume cleaned events already exist.
  --help          Show help message


Usage formula:
    .. code-block:: bash

        extract-nustar <OBSID> <SRC_REG> <BKG_REG> \
            [--indir DIR] [--outdir DIR] [--ra RA --dec DEC] [--no-pipe]

Usage examples:
    .. code-block:: bash

        extract-nustar 80902505002 source.reg bkg.reg
        extract-nustar 80902505002 source.reg bkg.reg --ra 178.231 --dec 15.572

Outputs:
    1. Same outputs as `nupipeline <https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/nupipeline.html>`_ (in ``./output/<OBSID>`` directory) and `nuproducts <https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/nuproducts.html>`_ (in ``./products/<OBSID>`` directory)
    2. Inside ``./products/<OBSID>`` directory, there will be ``FPMA/`` and ``FPMB/`` sub-directories with spectrum files that will be used in :py:class:`~xsnap.spectrum.SpectrumFit`, i.e. ``nu<OBSID>A01_sr.pha`` and ``nu<OBSID>B01_sr.pha``

Requirements:
  - `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ with `NuSTAR CALDB files <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_supported_missions.html>`_ installed