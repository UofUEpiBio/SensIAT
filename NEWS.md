# SensIAT (development version)

# SensIAT 0.3.0

* Conformed most names to standard of starting with verbs.
* Documentation fixes and updates.
* Removed unnecessarily exported functions.
* Updated single index model with norm 1 coefficients to allow for multiple iterations.
* Various bug fixes.
* Source code cleanup.


# SensIAT 0.2.0

* New options for fitting single index outcome models:
    + MAVE (Minimum Average Variance Estimation)
    + Fixed bandwidth, where the bandwidth is held constant and all coefficients are allowed to vary freely.

# SensIAT 0.1.1

* Bug fix for compiling on Debian and other systems without C++20.

# SensIAT 0.1.0

* Initial CRAN submission.
* Beta Release.
* Includes support for `'dnorm'` and `'K2_biweight'`(quartic) kernels.
