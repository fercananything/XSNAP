# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/)
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.23] - 2025-07-21
### Changed
- Edit the help docs for the scripts functions

## [0.1.22] - 2025-07-21
### Changed
- `SpectrumFit.get_fluxes()` now find every component that has nH in the model and set it to zero for unabsorbed flux calculations

## [0.1.21] - 2025-07-17
### Changed
- Bug fixes

## [0.1.20] - 2025-07-17
### Changed
- CSMAnalysis is now only mcmc-based fit

## [0.1.17-19] - 2025-07-17
### Changed
- Bug fixes

## [0.1.16] - 2025-07-15
### Changed
- `extract-nustar` now can be used with just parsing a background region, will automatically use the source region from the downloaded data

## [0.1.15] - 2025-07-15
### Changed
- Name of function in SpectrumManager, from `plot_flux_lc()` to `plot_flux()`, same with `plot_lumin_lc()`, `plot_phot_lc()`, and `plot_counts_lc()`

## [0.1.14] - 2025-07-15
### Added
- URLs in `pyproject.toml`

## [0.1.13] - 2025-07-15
### Added
- Upload to Github at [https://github.com/fercananything/XSNAP/](https://github.com/fercananything/XSNAP/)

## [0.1.12] - 2025-07-15
### Changed
- README.rst become README.md
- Wrote descriptions and installations guide in README.md
- Defined the extract-nustar script in pyproject.toml

## [0.1.11] - 2025-07-14
### Changed
- Edit SpectrumManager plot_counts_lc() to be grouped by instrument

## [0.1.10] – 2025-07-14
### Changed
- Edit SpectrumFit on simulate function error

## [0.1.9] – 2025-07-14
### Changed
- Edit SpectrumFit on simulate function to allow plotting the simulated data

## [0.1.8] – 2025-07-14
### Changed
- Edit SpectrumFit on simulate function for upper limits to use raw counts

## [0.1.7] – 2025-07-13
### Changed
- Fixed the fitting asymmetric freezing norm/exp to density mcmc error
- README.md become README.rst

## [0.1.6] – 2025-07-13
### Changed
- Edit the fitting asymmetric freezing norm/exp to density mcmc error

## [0.1.5] – 2025-07-13
### Changed
- Add the fitting asymmetric freezing norm/exp to density mcmc

## [0.1.4] – 2025-07-13
### Changed
- Add mute function to SpectrumFit by using xspec.Xset.chatter
- Add count rates for simulated spectrum
- Change fitting asymmetric error to be able to freeze the norm/exp

## [0.1.3] – 2025-07-13
### Changed
- SpectrumManager is inside spectrum.py now

## [0.1.2] – 2025-07-13
### Added
- Add NuSTAR extract script

### Changed
- extract-swift script is changed to be more flexible

## [0.1.1] – 2025-07-10
### Changed
- Incorporated distance error functions
- Add simulated function in SpectrumFit to count upper limit fluxes

## [0.1.0] – 2025-07-10
### Added
- First release: all packages, LICENSE, README.md, CHANGELOG.md, pyproject.toml
- CLI script for Chandra, Swift-XRT, and XMM extraction, make_region for source and background based off an evt file
- Basic XSPEC fitting and supernova CSM analysis.
- Published in test.pypi.org