.. _make-region::

**********************************************
Making Region Files (:code:`make-region`)
**********************************************

Pipeline module to generate DS9-compliant ICRS source and background regions.

This module can automate region creation, 
including contamination checking via `XIMAGE <https://heasarc.gsfc.nasa.gov/docs/xanadu/ximage/ximage.html>`_. 
Thus, ``make-region`` is a good preliminary for making source and background regions. 
If user has DS9 installed, ``make-region`` will open DS9 so further adjustments can be made.

Positional arguments:
  - ``evtfile``  :  Path to event FITS file.
  - ``ra, dec`` : Source coordinates in right ascension (RA) and declination in ICRS (decimal degrees).
  - ``r_in`` :    Inner radius of source region in arcsec (default is 2" for Chandra, 45" for NuSTAR, and 25" for others).
  - ``r_out`` :   Outer radius for background region in arcsec (default is 45" for Chandra, 125" for others).
  - ``outdir`` :  Output directory (default: same as ``evtfile`` directory).

Options:
  --ds9      Path to DS9 executable (if not available on ``which ds9``).
  --expimg   Path to exposure image file.
  --help     Show help message

Output files:
  1. ``source.reg`` : circular region at the specified (RA, dec) with radius ``r_in`` in ICRS.
  2. ``bkg.reg`` :  Annulus (if no contamination around source) with inner radius ``r_in`` and outer radius ``r_out`` or circular region with radius ``r_out`` in ICRS. ``make-region`` automatically relocates background if there is an `XIMAGE <https://heasarc.gsfc.nasa.gov/docs/xanadu/ximage/ximage.html>`_ detection to a cleaner spot within around 1' from the source
  3. ``source_physical.reg`` : source region in physical coordinates (will only be generated if DS9 is executable)
  4. ``bkg_physical.reg`` : background region in physical coordinates (will only be generated if DS9 is executable)

Usage formula:
    >>> make_region <evtfile> <ra> <dec> [r_in] [r_out] [outdir] [--ds9 DS9_PATH] [--expimg EXP_IMG]

Usage examples: 
    >>> make_region /path/to/event.fits 169.532 -32.456
    >>> make_region /path/to/event.fits --expimg /path/to/exposure.img 169.532 -32.456
    >>> make_region /path/to/event.fits --expimg /path/to/exposure.img 169.532 -32.456 25 150 ./output/

Optional requirements:
    - `DS9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`_ for more interactive region adjustments