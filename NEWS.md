# SensIAT (development version)

* Added generalized fitting support across marginal mean and within-group workflows, including flexible link/loss handling and GLM outcome-model support.
* Added new simulation and resampling capabilities, including intensity/outcome simulator factories and a parametric bootstrap workflow.
* Improved performance and numerical stability for generalized fitting and influence-term computation, including caching and additional integration-method support.
* Expanded tests, CI checks, and documentation, alongside multiple bug fixes discovered during broader generalized-model and bootstrap development.

# SensIAT 0.3.0

* Conformed most names to standard of starting with verbs.
* Documentation fixes and updates.
* Removed unnecessarily exported functions.
* Updated single index model with norm 1 coefficients to allow for multiple iterations.
* Various bug fixes.
* Source code cleanup.


# SensIAT 0.2.0

* New options for fitting single index outcome models:
    + `MAVE` (Minimum Average Variance Estimation)
    + Fixed bandwidth, where the bandwidth is held constant and all coefficients are allowed to vary freely.

# SensIAT 0.1.1

* Bug fix for compiling on Debian and other systems without C++20.

# SensIAT 0.1.0

* Initial CRAN submission.
* Beta Release.
* Includes support for `'dnorm'` and `'K2_biweight'`(quartic) kernels.
