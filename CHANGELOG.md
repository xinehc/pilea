# Changelog
## [1.2.4] - 2025-07-22
### Fixed
- Cap KMC threads to 128.


## [1.2.3] - 2025-07-11
### Changed
- Simplify database format to reduce computational time for highly redundant database.
- Temporary file `.sketch` -> `.kmc`.

## [1.2.2] - 2025-07-10
### Fixed
- Reduce peak memory usage for highly redundant database.


## [1.2.1] - 2025-06-30
### Fixed
- Add np.eps to R to prevent RuntimeWarning in rare cases.


## [1.2.0] - 2025-06-25
### Added
- Support file list as input for indexing.

### Changed
- Avoid reconstructing dict for kmer reassignment.

### Fixed
- Ignore errors for `rmtree`.


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
