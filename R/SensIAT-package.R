#' @keywords internal
"_PACKAGE"

#' @import dplyr

## usethis namespace: start
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.count
#' @importFrom generics prune
#' @importFrom ggplot2 autoplot
#' @importFrom glue glue
#' @importFrom KernSmooth dpill
#' @importFrom MASS glm.nb
#' @importFrom methods is
#' @importFrom methods new
#' @importFrom orthogonalsplinebasis evaluate
#' @importFrom orthogonalsplinebasis GramMatrix
#' @importFrom orthogonalsplinebasis integrate
#' @importFrom orthogonalsplinebasis OBasis
#' @importFrom orthogonalsplinebasis orthogonalize
#' @importFrom orthogonalsplinebasis SplineBasis
#' @importFrom pracma quadv
#' @importFrom purrr discard_at
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
#' @importFrom stats dbinom
#' @importFrom stats delete.response
#' @importFrom stats dnbinom
#' @importFrom stats dnorm
#' @importFrom stats dpois
#' @importFrom stats family
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats na.fail
#' @importFrom stats na.omit
#' @importFrom stats na.pass
#' @importFrom stats nlminb
#' @importFrom stats optim
#' @importFrom stats predict
#' @importFrom stats predict.glm
#' @importFrom stats qnbinom
#' @importFrom stats qnorm
#' @importFrom stats qpois
#' @importFrom stats rexp
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats terms
#' @importFrom stats time
#' @importFrom stats update.formula
#' @importFrom stats var
#' @importFrom survival coxph
#' @importFrom survival strata
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom tibble tibble
#' @importFrom tidyr complete
#' @importFrom utils head
#' @importFrom utils tail
#' @useDynLib SensIAT, .registration = TRUE
## usethis namespace: end
NULL

globalVariables(c('..id..', '..time..', '..outcome..', '..prev_outcome..',
                  '..visit_number..', '..delta_time..', '..prev_time..'))
