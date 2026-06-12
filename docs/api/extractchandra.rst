.. _extract-chandra:

**********************************************
Extract CXO Data (:code:`extract-chandra`)
**********************************************

Pipeline module to generate Chandra images, PSF maps, source/background
regions, and extracted spectra from a Chandra X-ray Observatory (CXO)
``evt2`` event file.

Automates the standard `CIAO <https://cxc.harvard.edu/ciao/>`_ workflow:

  1. Locate the user’s Conda environment.
  2. Locate the appropriate aspect solution (``.asol1.fits``) and mask
     (``.msk1.fits``) files.
  3. Convert the requested RA/Dec to SKY coordinates with ``dmcoords``
     when auto-generating regions.
  4. Run ``dmcopy`` to make a filtered image.
  5. Create a PSF map with ``mkpsfmap``.
  6. If source/background regions are not supplied, run ``wavdetect`` and
     auto-create physical/SKY regions.
  7. If source/background regions are supplied but are not physical regions,
     convert simple ``fk5``/``icrs`` circle or annulus regions to physical/SKY
     coordinates.
  8. Run ``specextract`` to build source and background spectra.

Auto-region behavior
====================

If ``--src-reg`` and ``--bkg-reg`` are not supplied:

  - ``--ra`` and ``--dec`` are required.
  - ``wavdetect`` is run on the image cutout.
  - If one or more ``wavdetect`` detections fall within
    ``--detect-match-radius`` of the input RA/Dec, the closest detection is
    used as the source center.
  - The selected center is converted to physical/SKY ``X,Y`` using
    ``dmcoords``.
  - A circular source region is written with radius ``--src-radius``.
  - A background annulus is written with inner radius ``--bkg-inner`` and
    outer radius ``--bkg-outer``.
  - By default, ``--bkg-inner`` is set to ``--src-radius``.
  - If a ``wavdetect`` detection falls inside the background annulus, the
    annulus is replaced by a circular background region with radius
    ``--bkg-circle-radius``.

Default auto-region values:

  - ``--detect-match-radius 5`` arcsec
  - ``--src-radius 2`` arcsec
  - ``--bkg-inner`` defaults to ``--src-radius``
  - ``--bkg-outer 45`` arcsec
  - ``--bkg-circle-radius`` defaults to ``--bkg-outer``

Existing-region behavior
========================

If ``--src-reg`` and ``--bkg-reg`` are supplied:

  - ``--ra`` and ``--dec`` are not required.
  - Existing physical regions are used directly.
  - Simple celestial ``fk5`` or ``icrs`` ``circle`` and ``annulus`` regions are
    converted to physical/SKY regions before running ``specextract``.
  - If no RA/Dec is supplied, the script makes a full-field image instead of a
    centered cutout.

Supported region shapes for automatic conversion:

  - ``circle(ra,dec,radius)``
  - ``annulus(ra,dec,r_in,r_out)``

Positional arguments
====================

``COND_ENV``
    Name of the Conda environment to check against ``$CONDA_DEFAULT_ENV``.

``EVT_FILE``
    Path to the Chandra ``evt2`` event file.

Options
=======

``--src-reg SRC_REG``
    Existing source region file.

``--bkg-reg BKG_REG``
    Existing background region file.

``--outdir OUTDIR``
    Output directory. Defaults to the event file directory.

``--ra RA``
    Source RA in degrees. Required only when auto-generating regions.

``--dec DEC``
    Source Dec in degrees. Required only when auto-generating regions.

``--src-radius SRC_RADIUS``
    Auto source radius in arcsec. Default: ``2``.

``--bkg-inner BKG_INNER``
    Auto background annulus inner radius in arcsec. Defaults to
    ``--src-radius``.

``--bkg-outer BKG_OUTER``
    Auto background annulus outer radius in arcsec. Default: ``45``.

``--bkg-circle-radius BKG_CIRCLE_RADIUS``
    Circular background radius in arcsec, used when a detection falls inside
    the proposed background annulus. Defaults to ``--bkg-outer``.

``--detect-match-radius DETECT_MATCH_RADIUS``
    Match radius in arcsec for choosing a ``wavdetect`` source as the region
    center. Default: ``5``.

``--wav-scales WAV_SCALES``
    ``wavdetect`` wavelet scales in image pixels. Default: ``"1 2 4 8"``.

``--wav-sigthresh WAV_SIGTHRESH``
    ``wavdetect`` significance threshold. Default: ``1e-6``.

``--cutout-size CUTOUT_SIZE``
    Image cutout size in SKY pixels when RA/Dec are supplied. Default:
    ``512``.

``--bin-size BIN_SIZE``
    Image bin size in SKY pixels. Default: ``1``.

``--emin EMIN``
    Minimum ACIS energy in eV. Default: ``500``.

``--emax EMAX``
    Maximum ACIS energy in eV. Default: ``8000``.

``--ccd-id CCD_ID``
    CCD ID filter. If omitted for ACIS energy-filtered data, defaults to
    ``7``.

``--no-energy-filter``
    Disable energy filtering. This should be used for HRC event files.

``--psf-energy PSF_ENERGY``
    Energy for ``mkpsfmap`` in keV. Default: ``1.4967``.

``--ecf ECF``
    Encircled counts fraction for ``mkpsfmap``. Default: ``0.90``.

Output files
============

Inside ``OUTDIR``:

  - ``evt_img.fits``
  - ``evt_psfmap.fits``
  - ``wavdetect_src.fits`` when auto-generating regions
  - ``source.reg``
  - ``bkg.reg``
  - ``spec{obsid}.arf``
  - ``spec{obsid}.rmf``
  - ``spec{obsid}.pi``
  - ``spec{obsid}.corr.arf``
  - ``spec{obsid}_bkg.arf``
  - ``spec{obsid}_bkg.rmf``
  - ``spec{obsid}_bkg.pi``
  - ``spec{obsid}_grp.pi``

.. note::

    Assuming that the event file is inside ``obsid/primary/event.fits``, the
    obsid is taken from the parent of the parent directory of the event file.

Usage formula
=============

Auto regions:

.. code-block:: bash

    extract-chandra <COND_ENV> <EVT_FILE> --ra <RA> --dec <DEC> [OPTIONS]

Existing regions:

.. code-block:: bash

    extract-chandra <COND_ENV> <EVT_FILE> --src-reg <SRC_REG> --bkg-reg <BKG_REG> [OPTIONS]

Usage examples
==============

HRC, auto regions:

.. code-block:: bash

    extract-chandra ciao-4.18 hrc_evt2.fits \
        --ra 10.83625 \
        --dec 41.30911 \
        --outdir . \
        --no-energy-filter \
        --cutout-size 512 \
        --bin-size 1

HRC, existing regions:

.. code-block:: bash

    extract-chandra ciao-4.18 hrc_evt2.fits \
        --src-reg source.reg \
        --bkg-reg bkg.reg \
        --outdir . \
        --no-energy-filter \
        --bin-size 1

ACIS, auto regions:

.. code-block:: bash

    extract-chandra ciao-4.18 acis_evt2.fits \
        --ra 10.83625 \
        --dec 41.30911 \
        --outdir . \
        --emin 500 \
        --emax 8000 \
        --cutout-size 512 \
        --bin-size 1

ACIS, existing regions:

.. code-block:: bash

    extract-chandra ciao-4.18 acis_evt2.fits \
        --src-reg source.reg \
        --bkg-reg bkg.reg \
        --outdir . \
        --emin 500 \
        --emax 8000 \
        --bin-size 1

Requirements:
  • `CIAO <https://cxc.harvard.edu/ciao/>`_ (including ``dmcopy``, ``mkpsfmap``, ``specextract``, ``wavdetect``, ``dmcoords``) installed.
  • ``CONDA_DEFAULT_ENV`` set to the requested environment name.
  • ``XSNAP`` must also be installed in ``CONDA_ENV``.
