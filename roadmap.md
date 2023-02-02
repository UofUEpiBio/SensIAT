# (WIP) Package development roadmap

To get the PCORI R package published on CRAN, we need to complete the following activities:

- Write the R functions.

- Document the package using `roxygen` or similar.

- Have `x` number of vignettes written with `rmarkdown` or `quarto` (not sure about the latter).

- Have a testing framework running with `tinytest` or `testthat` (I prefer the former.)

- Submit the R package and test it in `(Windows, OSX, Ubutnu) x (R current, R dev, R old rel)`.

Each one of those activities has complexities of their own. In my experience, setting up the testing framework could be the most time-consuming.

