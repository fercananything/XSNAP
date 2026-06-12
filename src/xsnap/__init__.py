"""
xray_pipeline
=============

Modular pipeline for extracting, stacking and fitting X-ray spectra
from Chandra, Swift-XRT and XMM-Newton.
"""
try:
    from .spectrum import SpectrumFit, SpectrumManager
    from .analysis import CSMAnalysis
    from .detect import SourceDetection
    from .temperature import TemperatureEstimator
except ImportError:
    SpectrumFit = None
    SpectrumManager = None
    CSMAnalysis = None
    SourceDetection = None
    TemperatureEstimator = None

__all__ = [ "SpectrumFit", "SpectrumManager", "CSMAnalysis",
            "SourceDetection", "TemperatureEstimator"
            ]
__version__ = "0.2.6"