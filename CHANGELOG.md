# Changelog
All notable changes to this project will be documented in this file.

## [0.2.9] - 2025-12-23
### Changed
- Update `extract-chandra` argument of `CCD_ID` not automatically `CCD_ID=7`. `CCD_ID` must be manually inputted


## [0.2.8] - 2025-12-23
### Changed
- Update `extract-swift` ARFFILE and BACKFILE name to not be a bug for `xsnap.SpectrumFit`

## [0.2.7] - 2025-12-07
### Changed
- Update `extract-chandra` to be more flexible with observations from ACIS / HRC.

## [0.2.6] - 2025-12-07
### Changed
- Description of `xsnap.TemperatureEstimator.compute_pl_fit()` as there has been some typos.

## [0.2.5] - 2025-07-30
### Changed
- Change contacts to support@xsnap.org and homepage is now in [https://xsnap.org](https://xsnap.org)

## [0.2.4] - 2025-07-30
### Changed
- Edit `__init__.py` to access the Python classes directly from XSNAP, e.g. `xsnap.SpectrumFit` instead of `xsnap.spectrum.SpectrumFit`

## [0.2.3] - 2025-07-30
### Changed
- `__init__.py` to access the Python classes directly from XSNAP, e.g. `xsnap.SpectrumFit` instead of `xsnap.spectrum.SpectrumFit`

## [0.2.2] - 2025-07-30
### Added
- First release in PyPI: all packages, LICENSE, README.md, CHANGELOG.md, pyproject.toml
- Python scripts for Chandra, Swift-XRT, and XMM extraction, make_region for source and background based off an event file
- Basic XSPEC fitting and supernova CSM analysis.
- Published in [test.pypi.org](https://test.pypi.org/p/xsnap) and [pypi.org](https://pypi.org/p/xsnap)
- Documentation is available [here](https://fercananything.github.io/XSNAP)