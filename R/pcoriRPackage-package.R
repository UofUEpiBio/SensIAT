#' @keywords internal
"_PACKAGE"

#' @import dplyr

## usethis namespace: start
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.count
#' @importFrom glue glue
#' @importFrom KernSmooth dpill
#' @importFrom MASS glm.nb
#' @importFrom methods is
#' @importFrom orthogonalsplinebasis evaluate
#' @importFrom orthogonalsplinebasis GramMatrix
#' @importFrom orthogonalsplinebasis integrate
#' @importFrom orthogonalsplinebasis OBasis
#' @importFrom orthogonalsplinebasis orthogonalize
#' @importFrom orthogonalsplinebasis SplineBasis
#' @importFrom pracma quadv
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom purrr pmap
#' @importFrom purrr reduce
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang :=
#' @importFrom rlang ensym
#' @importFrom rlang eval_tidy
#' @importFrom splines ns
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats dnbinom
#' @importFrom stats dnorm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats na.omit
#' @importFrom stats na.pass
#' @importFrom stats nlminb
#' @importFrom stats optim
#' @importFrom stats predict
#' @importFrom stats predict.glm
#' @importFrom stats rexp
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats terms
#' @importFrom stats time
#' @importFrom stats var
#' @importFrom survival coxph
#' @importFrom survival strata
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom tibble tibble
#' @importFrom tidyr complete
#' @importFrom utils head
#' @importFrom utils tail
#' @useDynLib pcoriRPackage, .registration = TRUE
## usethis namespace: end
NULL
