#' Fit the marginal mean model for generalize outcomes.
#'
#' This function uses adaptive Simpson's rule for integration and supports multiple
#' outcome model types (single-index, GLM, linear models) through the generic
#' `compute_SensIAT_expected_values` interface. The implementation includes
#' numerical stability improvements (exp(-μ) multiplication vs division) and
#' caching optimizations for repeated expected value computations.
#'
#' @inheritParams fit_SensIAT_marginal_mean_model
#' @param loss The loss function to use. Options are "lp_mse", "mean_mse", and "quasi-likelihood".
#' @param link The link function to use. Options are "identity", "log", and "logit".
#' @param term2_method Method for computing term2 influence components. Options are "fast" (default, optimized closure-based integrand) and "original" (standard implementation).
#'
#' @details
#' ## Integration Method
#' The function uses adaptive Simpson's quadrature (`pcoriaccel_integrate_simp`)
#' for numerical integration, which automatically adjusts step sizes based on
#' local function behavior.
#'
#' ## Outcome Model Compatibility
#' Unlike simulation code that assumes specific single-index model formulas, this
#' function supports any outcome model with a `compute_SensIAT_expected_values`
#' method, including:
#' - Single-index models (`fit_SensIAT_single_index_*_model`)
#' - Generalized linear models (GLM)
#' - Linear models (LM)
#' - Negative binomial models
#'
#' ## Performance Optimizations
#' - Expected values for term1 are cached (depend only on alpha, not beta)
#' - Weight functions use numerically stable exp(-μ) multiplication
#' - Fast method uses closure-based integration with reduced allocations
#'
#'
#' @examples
#' library(survival)
#'
#' data_with_lags <- SensIAT_example_data |>
#'     dplyr::group_by(Subject_ID) |>
#'     dplyr::mutate(
#'         ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
#'         ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
#'         ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
#'     )
#'
#' # Create the observation time intensity model
#' intensity.model <-
#'     coxph(Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
#'         data = data_with_lags |> dplyr::filter(.data$Time > 0)
#'     )
#'
#' # Create the observed outcome model
#' outcome.model <-
#'     fit_SensIAT_single_index_fixed_coef_model(
#'         Outcome ~ ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
#'         id = Subject_ID,
#'         data = data_with_lags |> filter(Time > 0)
#'     )
#' fit_SensIAT_marginal_mean_model_generalized(
#'     data = data_with_lags,
#'     time = data_with_lags$Time,
#'     id = data_with_lags$Subject_ID,
#'     alpha = 0,
#'     knots = c(60, 260, 460),
#'     outcome.model = outcome.model,
#'     intensity.model = intensity.model,
#'     loss = "lp_mse",
#'     link = "log",
#'     impute_data = \(t, df){
#'         data_wl <- df |>
#'             mutate(
#'                 ..prev_time.. = Time,
#'                 ..prev_outcome.. = Outcome,
#'                 ..delta_time.. = 0
#'             )
#'         extrapolate_from_last_observation(t, data_wl, "Time", slopes = c("..delta_time.." = 1))
#'     }
#' )
#' time <- data_with_lags$Time
#' id <- data_with_lags$Subject_ID
#'
#' @export
fit_SensIAT_marginal_mean_model_generalized <-
    function(data,
             time, #< 0-indexed time vector
             id, #< Integer vector of patient ids
             alpha,
             knots,
             outcome.model,
             intensity.model,
             impute_data,
             loss = c("lp_mse", "quasi-likelihood"),
             link = c("identity", "log", "logit"),
             spline.degree = 3L,
             ..., #< passed to impute_data
             BBsolve.control = list(
                 maxit = 1000,
                 tol = 1e-6
             ),
             term2_method = c("fast", "original")) {
        time <- enquo(time)
        time <- rlang::eval_tidy({{time}}, data)

        id <- enquo(id)

        Y <- dplyr::pull(data, rlang::f_lhs(formula(outcome.model)))


        # match arguments
        loss <- match.arg(loss)
        link <- match.arg(link)
        term2_method <- match.arg(term2_method)

        if (link == "identity") {
            return(
                fit_SensIAT_marginal_mean_model(
                    data = data,
                    id = !!id,
                    alpha = alpha,
                    knots = knots,
                    outcome.model = outcome.model,
                    intensity.model = intensity.model,
                    spline.degree = spline.degree,
                    ...
                )
            )
        }

        if (is(knots, "SplineBasis")) {
            base <- knots
        } else {
            knots <- c(
                rep(head(knots, 1), spline.degree),
                knots,
                rep(tail(knots, 1), spline.degree)
            )
            base <- SplineBasis(knots, order = spline.degree + 1L)
        }

        tmin <- base@knots[base@order]
        tmax <- base@knots[length(base@knots) - base@order + 1]

        if (loss == "lp_mse") {
            if (link == "log") {
                link.fun <- log
                inv.link <- exp
                d1.inv.link <- exp
                V <- GramMatrix(base)
                V.inv <- solve(V)

                W <- function(t, beta) {
                    B <- pcoriaccel_evaluate_basis(base, t)
                    mu <- as.numeric(B %*% beta)
                    as.vector((V.inv %*% B) * exp(-mu))
                }
            } else if (link == "logit") {
                link.fun <- function(mu) log(mu / (1 - mu))
                inv.link <- function(eta) exp(eta) / (1 + exp(eta))
                d1.inv.link <- function(eta) {
                    exp(eta) / ((1 + exp(eta))^2)
                }
                V <- GramMatrix(base)
                V.inv <- solve(V)

                W <- function(t, beta) {
                    B <- pcoriaccel_evaluate_basis(base, t)
                    eta <- as.numeric(B %*% beta)
                    # For logit link: dg/dz = 1/[z(1-z)]
                    # Evaluated at z = s(eta) = exp(eta)/(1+exp(eta))
                    # s(eta)(1-s(eta)) = exp(eta)/(1+exp(eta))^2
                    # So 1/[s(eta)(1-s(eta))] = (1+exp(eta))^2/exp(eta)
                    # W_1 = V^{-1} B(t) * (1+exp(eta))^2 / exp(eta)
                    # For numerical stability, rewrite as: V^{-1} B(t) * (1 + 2*exp(-eta) + exp(-2*eta))
                    as.vector((V.inv %*% B) * (1 + 2 * exp(-eta) + exp(-2 * eta)))
                }
            } else {
                stop("Unsupported link function for lp_mse loss.")
            }
            # } else if(loss == 'mean_mse'){
            #     if(link == 'log'){
            #         link.fun = log
            #         inv.link = exp
            #         d1.inv.link = exp
            #     } else if(link == 'logit'){
            #         link.fun = function(mu) log(mu / (1 - mu))
            #         inv.link = function(eta) exp(eta) / (1 + exp(eta))
            #         d1.inv.link = function(eta) {
            #             exp(eta) / ((1 + exp(eta))^2)
            #         }
            #     } else {
            #         stop("Unsupported link function for mean_mse loss.")
            #     }
        } else if (loss == "quasi-likelihood") {
            if (link == "identity") {
                link.fun <- function(mu) mu
                inv.link <- function(eta) eta
                d1.inv.link <- function(eta) rep(1, length(eta))
                V <- GramMatrix(base)
                V.inv <- solve(V)
                
                W <- function(t, beta) {
                    B <- pcoriaccel_evaluate_basis(base, t)
                    # For identity link: ds/dz = 1, so V_3 = V_1 (constant)
                    as.vector(V.inv %*% B)
                }
            } else if (link == "log") {
                link.fun <- log
                inv.link <- exp
                d1.inv.link <- exp
                
                # For log link, V_3 depends on beta, so we need numerical integration
                # V_3(beta) = int B(t) B(t)' * exp(B(t)'beta) dt
                # This requires numerical integration for each beta
                # W_3(t; beta) = V_3(beta)^{-1} B(t)
                
                # Placeholder V.inv (not used, but needed for function signature)
                V <- GramMatrix(base)
                V.inv <- solve(V)
                
                W <- function(t, beta) {
                    # Compute V_3(beta) via numerical integration
                    integrand <- function(s) {
                        Bs <- pcoriaccel_evaluate_basis(base, s)
                        eta_s <- as.numeric(crossprod(Bs, beta))
                        # Return B(s) B(s)' * exp(eta_s) as a vector (column-major)
                        as.vector(tcrossprod(Bs) * exp(eta_s))
                    }
                    
                    V3_result <- pcoriaccel_integrate_simp(integrand, tmin, tmax)
                    V3_vec <- V3_result$Q
                    n <- ncol(base)
                    V3 <- matrix(V3_vec, n, n)
                    V3.inv <- solve(V3)
                    
                    B <- pcoriaccel_evaluate_basis(base, t)
                    as.vector(V3.inv %*% B)
                }
            } else if (link == "logit") {
                link.fun <- function(mu) log(mu / (1 - mu))
                inv.link <- function(eta) exp(eta) / (1 + exp(eta))
                d1.inv.link <- function(eta) {
                    exp(eta) / ((1 + exp(eta))^2)
                }
                
                # For logit link: ds/dz = exp(z)/(1+exp(z))^2
                # V_3(beta) = int B(t) B(t)' * exp(B(t)'beta)/(1+exp(B(t)'beta))^2 dt
                
                # Placeholder V.inv (not used, but needed for function signature)
                V <- GramMatrix(base)
                V.inv <- solve(V)
                
                W <- function(t, beta) {
                    # Compute V_3(beta) via numerical integration
                    integrand <- function(s) {
                        Bs <- pcoriaccel_evaluate_basis(base, s)
                        eta_s <- as.numeric(crossprod(Bs, beta))
                        # ds/dz at z=eta_s
                        weight <- exp(eta_s) / (1 + exp(eta_s))^2
                        # Return B(s) B(s)' * weight as a vector (column-major)
                        as.vector(tcrossprod(Bs) * weight)
                    }
                    
                    V3_result <- pcoriaccel_integrate_simp(integrand, tmin, tmax)
                    V3_vec <- V3_result$Q
                    n <- ncol(base)
                    V3 <- matrix(V3_vec, n, n)
                    V3.inv <- solve(V3)
                    
                    B <- pcoriaccel_evaluate_basis(base, t)
                    as.vector(V3.inv %*% B)
                }
            } else {
                stop("Unsupported link function for quasi-likelihood loss.")
            }
        } else {
            stop("Unsupported loss function.")
        }

        id <- rlang::eval_tidy({{id}}, data)
        fu <- time > 0

        intensity <- estimate_baseline_intensity(
            intensity.model = intensity.model,
            data = data[fu, ]
        )
        exp_gamma <- predict(intensity.model, newdata = data[fu, ], type = "risk", reference = "zero")
        intensity_weights <- intensity$baseline_intensity * exp_gamma

        # Define the influence function.
        influence <- function(beta) {
            valid.time <- time[fu] >= tmin & time[fu] <= tmax
            weights <- purrr::map(time[fu][valid.time], W, beta = beta)

            # Cache expected values for term1 (they only depend on alpha, not beta)
            if (!exists(".expected_cache", envir = environment(influence))) {
                .expected_cache <- compute_SensIAT_expected_values(
                    model = outcome.model,
                    alpha = alpha,
                    new.data = data[fu, ][valid.time, ]
                )
                assign(".expected_cache", .expected_cache, envir = environment(influence))
            }
            expected <- get(".expected_cache", envir = environment(influence))

            term1.deviation.by.observation <-
                (Y[fu][valid.time] - expected$E_Yexp_alphaY / expected$E_exp_alphaY) / 
                intensity_weights[valid.time]

            term1.by.observation <- purrr::map2(weights, term1.deviation.by.observation, `*`)


            compute_term2_for_patient <- function(patient_id) {
                patient_data <- data[id == patient_id, ]
                if (term2_method == "fast") {
                    compute_term2_influence_fast(
                        patient_data = patient_data,
                        outcome_model = outcome.model,
                        base = base,
                        alpha = alpha,
                        marginal_beta = beta,
                        V_inv = V.inv,
                        tmin = tmin,
                        tmax = tmax,
                        impute_fn = impute_data,
                        inv_link = inv.link,
                        W = W
                    )
                } else {
                    compute_term2_influence_original(
                        patient_data = patient_data,
                        outcome_model = outcome.model,
                        base = base,
                        alpha = alpha,
                        marginal_beta = beta,
                        V_inv = V.inv,
                        tmin = tmin,
                        tmax = tmax,
                        impute_fn = impute_data,
                        inv_link = inv.link,
                        W = W
                    )
                }
            }

            term2.by.patient <- map(unique(id), compute_term2_for_patient)

            reduce(term2.by.patient, `+`) + reduce(term1.by.observation, `+`)
        }

        influence(rep(1 / ncol(base), ncol(base)))
        time <- system.time(
            solution <- BB::sane(
                par = rep(1 / ncol(base), ncol(base)),
                fn = influence,
                control = BBsolve.control
            )
        )
    }
