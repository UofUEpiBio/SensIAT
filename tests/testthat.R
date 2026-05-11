# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(SensIAT)

# Use the verbose debug reporter so CI shows every expectation and context.
# This helps identify exactly where tests are hanging or slowing down.
test_check("SensIAT", reporter = DebugReporter$new())
