# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(SensIAT)

# Suppress solver output from optimization routines in test output
# and use ProgressReporter for clear per-test feedback
options(
  # Suppress BB solver iteration messages
  BB.iter.warnings = FALSE,
  # Reduce clutter in CI output
  digits = 5
)

test_check("SensIAT", reporter = ProgressReporter$new(show_praise = TRUE))
