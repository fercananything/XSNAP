[![PyPI - Version](https://img.shields.io/pypi/v/xsnap)](https://pypi.org/project/xsnap/) 
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/xsnap)](https://pypi.org/project/xsnap/)
[![Downloads](https://img.shields.io/pepy/dt/xsnap)](https://pepy.tech/project/xsnap)

<p align="center">
  <picture>
    <img src="https://raw.githubusercontent.com/fercananything/xsnap/main/docs/_static/logo/xsnap_logo_icon_crop.jpeg" alt="My Project logo" width="200px" />
  </picture>
</p>

# XSNAP: X-ray Supernova Analysis Pipeline

XSNAP (X-ray Supernova Analysis Pipeline) is a Python-based pipeline module that automates every step of X-ray supernova data reduction and analysis, from raw event processing and region selection to spectral fitting. XSNAP provides dedicated standard data calibration and spectral extraction scripts for Chandra X-ray Observatory (CXO), Swift-XRT, XMM-Newton, and NuSTAR data.

XSNAP, with the help of PyXspec, is able to model and fit spectra using a wide range of astrophysical models (e.g., Thermal-Bremsstrahlung and Powerlaw). Additionally, XSNAP can generate photometric data through the fitted spectra. 

A follow-up analysis using the Thermal-Bremsstrahlung model can be made, specifically for Type II Supernova. From luminosity fitting to estimating Circumstellar Medium (CSM) densities and mass-loss rates of the supernova progenitors, XSNAP streamines the workflow so you can spend less time on rewriting each analysis manually.

Documentation of XSNAP can be seen in this [page](https://xsnap.org/)

More on analysis functions can be made upon requests :)

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

Currently, [XSNAP](https://pypi.org/p/xsnap) is available for download in the popular Python Package Index (PyPI).
```shell script
pip install xsnap
```

Additionally, XSNAP should be able to be downloaded by cloning this Github repository and run:
```shell script
git clone https://github.com/fercananything/XSNAP/
cd XNAP
python -m pip install .
```

## Dependencies
### Required Dependencies

XSNAP analysis depends heavily on two non-Python softwares:
* [HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html) - Version 6.35. Other recent versions should be compatible even if I have yet to test it. It is best to install all the packages; otherwise, install at least `NuSTAR`, `Swift`, all of `General-Use FTOOLS`, and all of `XANADU`.
* [HEASOFT's PyXspec](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/buildinstall.html) - Version 2.1.4 (or XSPEC - Version 12.15.0). Other recent versions should be compatible even if I have yet to test it. Additionally, PyXspec should be automatically installed when you install HEASOFT.

### Recommended Dependencies

While it's not necessarily required, it is strongly recommended to download these non-Python softwares:

* [Chandra Interactive Analysis of Observations (CIAO)](https://cxc.harvard.edu/ciao/download/index.html) - Version 4.17. CIAO is needed if you want to do the spectral extraction from CXO data. It is recommended to install CIAO using the `conda create` command, i.e. install on a different Python/Conda virtual environment. This is to seperate HEASOFT (and XSPEC) with CIAO and avoid clashes between modules. 
* [XMM Science Analysis System (SAS)](https://www.cosmos.esa.int/web/xmm-newton/sas-download) - Version 22.1. However, other recent versions should still be compatible. SAS is needed if you want to do data calibration and spectral extraction for XMM-Newton. A few extra steps for SAS installation can be found [here](https://www.cosmos.esa.int/web/xmm-newton/sas-thread-startup#).
* [HEASARC Calibration Database (CALDB)](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/install.html) - Version 2009 Aug 04. The HEASARC CALDB is needed if you want to do data calibration and spectral extraction for Swift-XRT and NuSTAR.
* [CALDB Files for Swift-XRT and NuSTAR](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_supported_missions.html). In addition to the CALDB, the CALDB files are needed to be downloaded too. These files are needed if you want to do data calibration and spectral extraction for Swift-XRT and NuSTAR.

_Keep in mind, without these softwares, you are only able to import the spectra fitting and analysis modules. These softwares help with the scripts dealing for data calibration and spectral extraction._

### Optional Dependencies

This software is completely optional and has minimal impact on the user experience.
* [DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9) - Version 4.1 and above. DS9 is needed to help user's interactivity in making region files.

## How to use the module

XSNAP is organized into two main parts: command-line scripts (where users can invoke on the shell or jupyter notebook) and a built-in module or Python API (where you can import functions and classes).

There are six scripts available for users to run:
| Script             | Description                                              |
|--------------------|----------------------------------------------------------|
| `extract-chandra`  | Calibrate & extract spectrum from Chandra observations. |
| `extract-swift`    | Calibrate & extract spectrum from Swift-XRT (PC/WT mode available).     |
| `swift-stack-pc`   | Bin & stack Swift-XRT PC-mode data (default 1-day bins). |
| `extract-xmm`      | Calibrate & extract spectrum from XMM-Newton.           |
| `extract-nustar`      | Calibrate & extract spectrum from NuSTAR.           |
| `make-region`      | Generate ICRS source/background region files. (Physical region files will also be made if user has DS9)       |

A short tutorial on how to use XSNAP is available in jupyter notebooks [here](https://github.com/fercananything/XSNAP/tree/main/notebook)

## Problems and Questions
If you encounter a bug, or would like to make a feature request, please use the GitHub
[issues](https://github.com/fercananything/XSNAP/issues) page.

In addition, if you have further questions, feel free to send me an email at support@xsnap.org