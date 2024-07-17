#' PCORI Example Data
#'
#' A simulated dataset for use in the PCORI tutorial, testing and documentation.
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
"PCORI_example_data"
