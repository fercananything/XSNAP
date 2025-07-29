.. _spectrum:

*************************************
Spectrum (:code:`xsnap.spectrum`)
*************************************

The spectrum objects in XSNAP are designed to assist with fitting, modeling, and 
analyzing x-ray spectra from observational data files. XSNAP provides two main classes:

- :py:class:`~xsnap.spectrum.SpectrumFit` is intended for fitting and modeling individual spectra. It allows users to extract physical parameters such as fluxes and luminosities from a single spectrum through robust fitting with `PyXspec <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/index.html>`_.

- :py:class:`~xsnap.spectrum.SpectrumManager` is designed for handling a collection of spectra. It streamlines the analysis of multiple spectra, enabling users to efficiently plot light curves, track the evolution of spectral parameters over time, and perform bulk operations on datasets.


.. toctree::
   :maxdepth: 1
   :hidden:

   xsnap.spectrum.SpectrumFit <spectrum.SpectrumFit>
   xsnap.spectrum.SpectrumManager <spectrum.SpectrumManager>