# (WIP) Package development roadmap

To get the PCORI R package published on CRAN, we need to complete the following activities:

- Write the R functions and document the package using `roxygen` or similar.

- Have a testing framework running with [`tinytest`](https://github.com/markvanderloo/tinytest) or [`testthat`](https://testthat.r-lib.org/) (I prefer the former.)

- Have `x` number of vignettes written with `rmarkdown` or `quarto` (not sure about the latter).

- Submit the R package and test it in `(Windows, OSX, Ubutnu) x (R current, R dev, R old rel)`.

Each one of those activities has complexities of their own. In my experience, setting up the testing framework could be the most time-consuming. Let's think about each one in detailed:

## Write the R functions and document them

Since there are two general methods, I see two sets of functions:

1. The main model fit function, and

2. The methods associated with them, in particular: `print()`, `plot()`, `summary()`, `coef()`, `predict()`, and `vcov()`.

The methods for the resulting outputs provide a great deal of flexibility as other R functions can grab those to compute other things. One example is the `confit()` function.

I estimate that it may take about two months to write the two main fit functions and one additional month to write the S3 methods. Considering that, in the ideal world, writing functions should go along with documentation. I would incorporate an additional month for documentation. Furthermore, CRAN is requesting in-depth documentation, including:

- A short description.

- Details: Talk about the method and layout formulas (if possible).

- Returns: In-depth description of the output.

- Examples: having one or more examples of the function.

## Preparing the tests

Testing is fundamental, and it can take a lot of work. I would usually create not one but a handful of tests evaluating (i) typical cases, (ii) borderline cases, and (iii) potential errors (catching errors in data). Building tests requires generating artificial data for which we know exactly the expected results, fit the model, and hope that the model behaves as expected. Because of this complexity, I would allocate two months of work to this activity.


## Write `x` vignettes

Vignettes are of great help in R packages. Although not as complicated as a paper, the perfect vignette will have an extended example from the start to the end of the package. For instance, loading data, preparing the data (if needed), fitting the model, and doing the post-estimation analyses. I would allocate one month of work for each vignette.

## CRAN submission

If we followed all the previous steps to the detail, submitting an R package should be straightforward. At this point, we should ensure the package runs on multiple OS and versions of R (at least R-dev and R-release). Being conservative, it should take about one month to get this bit done (thinking about two to three submissions, just like in the paper review process.)

# Final count

- Write the R functions and document them: 3 months

- Preparing the tests: 2 months.

- Write `x` vignettes: Assuming 2, 2 months.

- CRAN submission: 1 month.

Total, 8 months. Providing some wiggle room (and accounting for the usual under estimation of time,) I would say 8 x 1.5 = 12 months (magic, but unexpected).
