# XSNAP: X-ray Supernova Analysis Pipeline

XSNAP (X-ray Supernova Analysis Pipeline) is a Python-based pipeline module that automates every step of X-ray supernova data reduction and analysis, from raw event processing and region selection to spectral fitting. XSNAP provides dedicated standard data calibration and spectral extraction scripts for Chandra X-ray Observatory (CXO), Swift-XRT, XMM-Newton, and NuSTAR data.

XSNAP, with the help of PyXspec, is able to model and fit spectra using a wide range of astrophysical models (e.g., Thermal-Bremsstrahlung and Powerlaw). Additionally, XSNAP can generate photometric data through the fitted spectra. 

A follow-up analysis using the Thermal-Bremsstrahlung model can be made, specifically for Type II Supernova. From luminosity fitting to estimating Circumstellar Medium (CSM) densities and mass-loss rates of the supernova progenitors, XSNAP streamines the workflow so you can spend less time on rewriting each analysis manually.

More analysis functions can be made upon requests :)

## Contents

1. [Introduction](#xsnap-x-ray-supernova-analysis-pipeline)  
2. [Installation](#installing-xsnap)  
3. [Dependencies](#required-dependencies)  
4. [Usage and Example](#how-to-use-the-module)  
   - [CLI Scripts](#command-line-scripts)  
   - [Python API](#built-in-module--python-api)  
5. [Problems and Questions](#problems-and-questions) 

## Installing XSNAP

We strongly recommend that you make use of Python virtual environments, or (even better) Conda virtual environments when installing XSNAP. 

Currently, XSNAP is in its development phase. It is available for download at the Testing Python Package Index (TestPyPI).
```shell script
pip install -i https://test.pypi.org/simple/ xsnap
```

<!-- 
XSNAP is available on the popular Python Package Index (PyPI), and can be installed like this:
```shell script
pip install xsnap
```
-->

## Required Dependencies

XSNAP analysis depends heavily on two non-Python softwares:
* [HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html) - Version 6.35. Other recent versions should be compatible even if I have yet to test it.
* [HEASOFT's PyXspec](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/buildinstall.html) - Version 2.1.4 (or XSPEC - Version 12.15.0). Other recent versions should be compatible even if I have yet to test it. Additionally, PyXspec should be automatically installed when you install HEASOFT.


## Recommended Dependencies

While it's not required, it is recommended to download these non-Python softwares:

* [Chandra Interactive Analysis of Observations (CIAO)](https://cxc.harvard.edu/ciao/download/index.html) - Version 4.17. CIAO is needed if you to do the spectral extraction from CXO data. It is recommended to install CIAO using the `conda create` command, i.e. install on a different Python/Conda virtual environment. This is to seperate HEASOFT (and XSPEC) with CIAO and avoid clashes between modules. 
* [HEASARC Calibration Database (CALDB)](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/install.html) - Version: 2009 Aug 04. The HEASARC CALDB is needed if you want to do data calibration and spectral extraction for Swift-XRT and NuSTAR.
* [CALDB Files for Swift-XRT and NuSTAR](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_supported_missions.html). In addition to the CALDB, the CALDB files are needed to be downloaded too. These files are needed if you want to do data calibration and spectral extraction for Swift-XRT and NuSTAR.

_Keep in mind, without these softwares, you are only able to import the spectra fitting and analysis modules. These softwares help with the scripts dealing for data calibration and spectral extraction._

## Optional Dependencies

This software is completely optional and has minimal impact on the user experience.
* [DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9) - Version 4.1 and above. DS9 is needed to help user's interactivity in making region files.

## How to use the module

XSNAP is organized into two main parts: command-line scripts (where users can invoke on the shell or jupyter notebook) and a built-in module or Python API (where you can import functions and classes).

### Command-line scripts

There are five scripts available for users to run:
| Script             | Description                                              |
|--------------------|----------------------------------------------------------|
| `extract-chandra`  | Calibrate & extract spectrum from Chandra observations. |
| `extract-swift`    | Calibrate & extract spectrum from Swift-XRT (PC/WT mode available).     |
| `swift-stack-pc`   | Bin & stack Swift-XRT PC-mode data (default 1-day bins). |
| `extract-xmm`      | Calibrate & extract spectrum from XMM-Newton.           |
| `extract-nustar`      | Calibrate & extract spectrum from NuSTAR.           |
| `make-region`      | Generate ICRS source/background region files. (Physical region files will also be made if user has DS9)       |

For details on how to use these scripts, you can try the `-help` argument.

For example:
```shell script
$ make-region -help

usage: make-region [-h] [--ds9 DS9] [--expimg EXPIMG]
                   evtfile ra dec [r_in] [r_out] [outdir]

Pipeline module to generate DS9-compliant ICRS source and background regions.

This script produces two region files for automated processing:

  1. source.reg
     • Circular region at the specified RA/Dec with radius r_in (″).
  2. bkg.reg
     • First-pass annulus from r_in to r_out (″).
     • If an XIMAGE detection (excluding the source) overlaps, relocates the
       background to a single circle of radius r_out (″) at the first clean spot
       ≥1′ from the source.

Command-line usage:
    make_region <evtfile> <ra> <dec> [r_in] [r_out] [outdir]
                [--ds9 DS9_PATH] [--expimg EXP_IMG]

Positional arguments:
  evtfile    Path to event FITS file.
  ra, dec    Source coordinates in ICRS (decimal degrees).
  r_in       Inner radius of source region (arcsec; defaults telescope-specific).
  r_out      Outer radius for background region (arcsec; defaults telescope-specific).
  outdir     Output directory (default: same as evtfile’s directory).

Options:
  --ds9      Path to DS9 executable (if not on $PATH or in $DS9/_PATH).
  --expimg   Exposure-map FITS file for contamination checks.
```

Then, you can run `make-region`, for example,
```shell script
make-region event.evt --expimg exposure_image.img 158 +28
```

### Built-in module / Python API

There are multiple functions and classes available. Below is an example on how to fit a spectrum.

```Python
# 1. Import xsnap
import xsnap
explosion_mjd = 60442.40 # Supernova explosion time in MJD

# 2. Extract spectrum and fit
spec = spectrum.SpectrumFit() 

# Load spectrum
spec.load_data('grp.pha')

# The same as rebin in XSPEC
spec.set_rebin(5, 1)

# Ignore energies outside 0.3 - 10.0 keV
spec.ignore("**-0.3 10.0-**")

# Set XSPEC model
spec.set_model( 
    "tbabs*ztbabs*bremss", # model args
    zTBabs_nH=2.7, 
    zTBabs_Redshift=0,
    bremss_kT=32,
    bremss_norm=0
) 
# TBabs.nH can be automatically defined by module using HEASOFT and the RA_OBJ and DEC_OBJ from spectrum file header

# Fit the spectrum and plot it (will have xAxis in keV as default)
spec.fit()
spec.set_plot("ldata", device="/svg")

# Fit the parameters and get their uncertainty
spec.get_params("1.0 2 5")

# Generate the energy and photon flux from the fitted models
df_flux = spec.get_fluxes()

# Get observation time of the spectrum
spec.get_time() 

# Get the count rate of the spectrum
spec.get_counts()

# Generate luminosity based on distance in Mpc (also available in redshift)
spec.get_lumin(df_flux['unabsorbed'], distance=14.1) 

# 3. Plot light curve
manager = spectrum.SpectrumManager(tExplosion=explosion_mjd)

# Load the SpectrumFit class to the SpectrumManager class
manager.load(spec)

# Generate a flux plot light curve in matplotlib
manager.plot_flux_lc()

# Generate a luminosity plot light curve in matplotlib
manager.plot_lumin_lc()
```

A notebook example will be made available in the near future.


## Problems and Questions
If you encounter a bug, or would like to make a feature request, or any questions, feel free to send me an email at ferdinand.1238073@gmail.com

<!-- If you encounter a bug, or would like to make a feature request, please use the GitHub
[issues](https://github.com/DavidT3/XGA/issues) page, it really helps to keep track of everything.

However, if you have further questions, or just want to make doubly sure I notice the issue, feel free to send
me an email at turne540@msu.edu -->