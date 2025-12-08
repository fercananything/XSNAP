"""
xray_pipeline
=============

Modular pipeline for extracting, stacking and fitting X-ray spectra
from Chandra, Swift-XRT and XMM-Newton.
"""
from .spectrum import SpectrumFit, SpectrumManager
from .analysis import CSMAnalysis
from .detect import SourceDetection
from .temperature import TemperatureEstimator

__all__ = [ "SpectrumFit", "SpectrumManager", "CSMAnalysis",
            "SourceDetection", "TemperatureEstimator"
            ]
__version__ = "0.2.6"