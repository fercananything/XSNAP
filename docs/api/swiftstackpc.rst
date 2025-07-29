.. _swift-stack-pc::

###########################################################
Stacking Swift/XRT Data (PC Mode) (:code:`swift-stack-pc`)
###########################################################

Pipeline module to run Swift/XRT processing in PC Mode across multiple epochs and stack time-binned products.

Automates the stacking workflow:
  1. Executes `xrtpipeline <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/help/xrtpipeline.html>`_ for each epoch listed in a text file, producing calibrated event and image files.
  2. Bins the resulting event (``.evt``) and image (``.img``) files into time intervals (days since explosion) defined by user.
  3. Stacks events and images in each bin into dedicated subdirectories under ``./stacks/``, named by their first–last epochs/observation id.
  4. Updates the stacked event files ``MJD-BEG`` and ``MJD-END`` headers.

.. caution::
    The user **must make sure** to run this script on the parent of ``"sources"`` directory and 
    the ``"sources"`` directory has the downloaded data with ``ObsId`` as its child directory. E.g.

    .. code-block:: bash

        parent
        ├── sources
        │   ├── ObsId1 # downloaded from Swift/XRT raw data
        │   ├── ObsId2
        |   ├── ObsId3
        |   └── ObsId4


    When ``swift-stack-pc`` is run inside parent directory, it will result in

    .. code-block:: bash

        parent
        ├── sources
        │   ├── ObsId1
        │   ├── ObsId2
        |   ├── ObsId3
        |   └── ObsId4
        ├── output
        │   ├── ObsId1
        │   ├── ObsId2
        |   ├── ObsId3
        |   └── ObsId4
        ├── stack
        │   ├── ObsId1_to_ObsId2
        │   └── ObsId3_to_ObsId4


Arguments:
  --epochs         File listing one epoch/obsid identifier per line.
  --bin-size       Time-bin width in days (default: ``1.0``).
  --explosion-mjd  Supernova time of explosion in MJD.
  --src-ra         Source ICRS RA (decimal degrees).
  --src-dec        Source ICRS Dec (decimal degrees).
  --pn-tra         (Optional) PN attitude RA override (default: ``POINT``).
  --pn-dec         (Optional) PN attitude Dec override (default: ``POINT``).
  --help           Show help message

Output:
   - The same output of `xrtpipeline <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/help/xrtpipeline.html>`_ in ``./output/`` directory.
   - ``./stacks/`` directory with their first–last epochs/observation id to be the subdirectories.
   - Inside the subdirectories of ``./stacks/``, there will be a cleaned evt file (format: `evt_{obsid1}to{obsid2}.evt`) and exposure image file (format: `img_{obsid1}to{obsid2}.img`).

Usage formula:
    .. code-block:: bash

        swift-stack-pc [--epochs] [--bin-size] [--explosion-mjd] \
            [--src-ra] [--src-dec] [--pn-tra] [--pn-dec]

Usage examples:
    .. code-block:: bash

        swift-stack-pc --epochs source_epochs.txt \
            --bin-size 2.0 \
            --explosion-mjd 59000.5 \
            --src-ra <RA> \
            --src-dec <Dec> \

        swift-stack-pc --epochs ./path/to/epochs.txt \
            --explosion-mjd 59000.5 \
            --src-ra <RA> \
            --src-dec <Dec> \
            

    

Requirements:
  - `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ with `Swift/XRT CALDB files <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_supported_missions.html>`_ installed