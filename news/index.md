# Changelog

## SensIAT (development version)

## SensIAT 0.3.0

CRAN release: 2025-09-05

- Conformed most names to standard of starting with verbs.
- Documentation fixes and updates.
- Removed unnecessarily exported functions.
- Updated single index model with norm 1 coefficients to allow for
  multiple iterations.
- Various bug fixes.
- Source code cleanup.

## SensIAT 0.2.0

CRAN release: 2025-08-18

- New options for fitting single index outcome models:
  - `MAVE` (Minimum Average Variance Estimation)
  - Fixed bandwidth, where the bandwidth is held constant and all
    coefficients are allowed to vary freely.

## SensIAT 0.1.1

CRAN release: 2024-11-17

- Bug fix for compiling on Debian and other systems without C++20.

## SensIAT 0.1.0

CRAN release: 2024-11-01

- Initial CRAN submission.
- Beta Release.
- Includes support for `'dnorm'` and `'K2_biweight'`(quartic) kernels.
