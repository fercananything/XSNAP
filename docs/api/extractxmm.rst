.. _extract-xmm::

**********************************************
Extract XMM-Newton Data (:code:`extract-xmm`)
**********************************************

Pipeline module to extract and group XMM-Newton EPIC (PN, MOS1, MOS2) spectra.

Automates the standard SAS (Science Analysis System) workflow:
  1. Build ``CCF``, ingest the ``ODF``, and set ``SAS`` environment variables.
  2. Locate the ``PPS`` directory alongside the ``ODF`` and prepare event files.
  3. Extract source and background spectra.
  4. Generate response files (RMF/ARF) and group spectra for PN, MOS1, and MOS2.
  5. Move final grouped spectra to the specified output directory.

.. caution::
    The user running this script **must make sure** that ``PPS`` and ``ODF`` directory have the same parent directory, e.g.,

    .. code-block:: bash

        OBSID
        ├── PPS/
        └── ODF/


Positional arguments:
  - ``ODF_DIR``    :    Path to the XMM-Newton ODF directory.
  - ``source.reg``  :   Source region in physical coordinates (filename or literal physical “circle(x,y,r)”).

    .. note::
        
        Keep in mind that the region files or literal can only be in circle(x, y, r) or annulus(x, y, r_in, r_out). 
        Also for literal, must be in string or "quoted".
  - ``PNbkg.reg``    :  PN background region in physical coordinates (filename or literal).
  - ``MOS1bkg.reg``   : Optional MOS1 background region in physical coordinates (defaults to ``PNbkg.reg`` if omitted).
  - ``MOS2bkg.reg``  :  Optional MOS2 background region in physical coordinates (defaults to ``PNbkg.reg`` if omitted).
  - ``OUTDIR``       :  Destination for grouped spectra in physical coordinates (defaul to ``PPS`` directory).

Options:
    --help      Show help message

Usage formula:
    >>> extract-xmm <ODF_DIR> <source.reg> <PNbkg.reg> [MOS1bkg.reg] [MOS2bkg.reg] [OUTDIR]

Usage examples:
    >>> extract-xmm ./0882480901/odf/ sou_PN_physical.reg bkg_PN_physical.reg \ 
    ...        bkg_MOS1_physical.reg bkg_MOS2_physical.reg ./0882480901
    >>> extract-xmm ./0882480901/odf/ sou_physical.reg bkg_physical.reg

Requirements:
  - `SAS (Science Analysis System) <https://www.cosmos.esa.int/web/xmm-newton/what-is-sas>`_ must be installed.
  - ``PPS`` and ``ODF`` directory must have the same parent directory
