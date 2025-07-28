

What is XSNAP?
===================

XSNAP (X-ray Supernova Analysis Pipeline) is a Python-based pipeline module that automates every step of X-ray supernova data reduction and analysis, from raw event processing and region selection to spectral fitting. XSNAP provides dedicated standard data calibration and spectral extraction scripts for Chandra X-ray Observatory (CXO), Swift-XRT, XMM-Newton, and NuSTAR data.

XSNAP, with the help of PyXspec, is able to model and fit spectra using a wide range of astrophysical models (e.g., Thermal-Bremsstrahlung and Powerlaw). Additionally, XSNAP can generate photometric data through the fitted spectra. 

A follow-up analysis using the Thermal-Bremsstrahlung model can be made, specifically for Type II Supernova. From luminosity fitting to estimating Circumstellar Medium (CSM) densities and mass-loss rates of the supernova progenitors, XSNAP streamines the workflow so you can spend less time on rewriting each analysis manually.

More analysis functions can be made upon requests (view `Development <development>`_ for more details)


