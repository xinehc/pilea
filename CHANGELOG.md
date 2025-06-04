# Changelog
## [1.1.2] - 2025-06-04
### Changed
- Rename parameter `--min-dept` to `--min-cove` and column `depth` to `coverage`.

### Fixed
- Fix docstring.


## [1.1.1] - 2025-05-31
### Fixed
- Fix a bug causing error in database merging (`-d`).


## [1.1.0] - 2025-05-25
### Added
- Add database version control.


## [1.0.0] - 2025-05-20
### Added
- Add alphabet checker for indexing.

### Changed
- Simplify parameter `--min-cont`.
- Reduce time and memory requirement for indexing.


## [0.2.3] - 2025-05-11
### Fixed
- Fix non-existence directory for `wget`.
- Fix file/sample mismatch caused by sorting.
- Fix Zenodo record id.


## [0.2.2] - 2025-04-18
### Changed
- Remove `.txt` for sketches.


## [0.2.1] - 2025-04-15
### Changed
- Refine criteria for existence checking of genomes.
- Simplify database format.


## [0.2.0] - 2025-04-08
### Added
- Add idr filter and loess-based qc bias correction.
### Changed
- Use non-duplicated sketches for fitting.


## [0.1.0] - 2025-03-11
### Added
- First release.
