#' SensIAT Example Data
#'
#' A simulated dataset for use in the SensIAT tutorial, testing and documentation.
#'
#' @format A data frame with 779 rows and 4 variables consisting of 200
#' simulated patients.  Each row in the data represents a visit for the patient.
#' The columns are:
#' \describe{
#'      \item{Subject_ID}{A unique identifier for each patient.}
#'      \item{Visit}{The ordinal number of the visit for the patient.  Baseline observation is 0.}
#'      \item{Time}{The time of the visit in days, since baseline.}
#'      \item{Outcome}{The outcome of interest.}
#' }
"SensIAT_example_data"


#' @describeIn SensIAT_example_data A simulated dataset with both treatment and control groups.
#' @format A data frame with 779 rows and 4 variables consisting of 200
#' simulated patients.  Each row in the data represents a visit for the patient.
#' The columns are:
#' \describe{
#'      \item{Subject_ID}{A unique identifier for each patient.}
#'      \item{Visit}{The ordinal number of the visit for the patient.  Baseline observation is 0.}
#'      \item{Time}{The time of the visit in days, since baseline.}
#'      \item{Outcome}{The outcome of interest.}
#'      \item{Treatment_group}{Treatment or control group.}
#' }
"SensIAT_example_fulldata"
