# Changelog

All notable changes to the CNMA package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned Features
- Meta-regression implementation
- Bayesian analysis backends (JAGS, Stan)
- Transportability weighting
- Publication bias diagnostics
- Interactive visualizations
- HTML report generation

## [1.0.0] - 2025-01-XX

### Added
- Initial release of CNMA package
- Core frequentist network meta-analysis via netmeta
- `run_cnma_analysis()` main analysis function
- `cnma_quickstart()` convenience function
- `setup_cnma()` comprehensive configuration system
- `validate_cnma_input()` input validation
- `cnma_clean_data()` data cleaning utilities
- `simulate_cnma_data()` demo data generation
- S3 methods: `print.cnma()`, `summary.cnma()`, `print.cnma_config()`
- Parallel processing support with future framework
- Caching mechanism for repeated analyses
- Comprehensive test suite (90+ tests)
- Professional documentation and README
- DESCRIPTION, NAMESPACE, and package infrastructure

### Changed
- Refactored monolithic 701-line file into 8 modular components
- Improved error messages with helpful hints
- Enhanced cache key generation with timestamp and random characters
- Better input validation with edge case handling

### Infrastructure
- Modular R/ directory structure:
  - `R/aaa-package.R` - Package documentation
  - `R/utilities.R` - Utility functions
  - `R/parallel.R` - Parallel processing
  - `R/config.R` - Configuration system
  - `R/validation.R` - Data validation
  - `R/data-simulation.R` - Data generation
  - `R/analysis.R` - Main analysis
  - `R/methods.R` - S3 methods
- Comprehensive test suite in `tests/testthat/`
- Professional documentation (README, NEWS, CONTRIBUTING)
- Git ignore rules and build configuration

### Documentation
- Comprehensive README with examples and usage guide
- Roxygen2 documentation for all exported functions
- CONTRIBUTING.md with development guidelines
- NEWS.md for tracking changes
- Clear roadmap for future development

### Known Limitations
- Bayesian analysis not yet implemented (planned for v1.2)
- Meta-regression not yet implemented (planned for v1.1)
- Transportability weighting not yet implemented (planned for v2.0)
- Advanced diagnostics not yet implemented (planned for v1.1-1.2)
- Interactive visualizations not yet implemented (planned for v1.2)

## Version Numbering

- **Major version** (X.0.0): Breaking changes to public API
- **Minor version** (1.X.0): New features, backward compatible
- **Patch version** (1.0.X): Bug fixes, backward compatible

---

For more details on each release, see [NEWS.md](NEWS.md).

[Unreleased]: https://github.com/mahmood726-cyber/HFN786/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/mahmood726-cyber/HFN786/releases/tag/v1.0.0
