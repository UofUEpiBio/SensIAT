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
             term2_method = c("fast", "original"),
             use_expected_cache = TRUE) {
        time.var <- enquo(time)
        time <- rlang::eval_tidy({{time.var}}, data)

        id <- enquo(id)

        Y <- dplyr::pull(data, rlang::f_lhs(formula(outcome.model)))


        # match arguments
        loss <- match.arg(loss)
        link <- match.arg(link)
        term2_method <- match.arg(term2_method)

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
            if (link == "identity") {
                # link.fun <- function(mu) mu
                inv.link <- function(eta) eta
                # d1.inv.link <- function(eta) rep(1, length(eta))
                V <- GramMatrix(base)
                V.inv <- solve(V)

                W <- function(t, beta) {
                    B <- pcoriaccel_evaluate_basis(base, t)
                    # For identity link: ds/dz = 1, so weight function is constant
                    as.vector(V.inv %*% B)
                }
            } else if (link == "log") {
                # link.fun <- log
                inv.link <- exp
                # d1.inv.link <- exp
                V <- GramMatrix(base)
                V.inv <- solve(V)

                W <- function(t, beta) {
                    B <- pcoriaccel_evaluate_basis(base, t)
                    mu <- as.numeric(B %*% beta)
                    as.vector((V.inv %*% B) * exp(-mu))
                }
            } else if (link == "logit") {
                # link.fun <- function(mu) log(mu / (1 - mu))
                inv.link <- function(eta) exp(eta) / (1 + exp(eta))
                # d1.inv.link <- function(eta) {exp(eta) / ((1 + exp(eta))^2)}
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
                    as.vector((V.inv %*% B) * (exp(eta) + 2 + exp(-eta)))
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
                # link.fun <- function(mu) mu
                inv.link <- function(eta) eta
                # d1.inv.link <- function(eta) rep(1, length(eta))
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
        unique_ids <- unique(id)

        # tmin assumed to be > 0 so that all included.obs imply followup.
        if (tmin == 0) {
            rlang::warn("tmin must be > 0")
        }
        included.obs <- time >= tmin & time <= tmax

        intensity <- estimate_baseline_intensity(
            intensity.model = intensity.model,
            data = data[included.obs, ]
        )
        exp_gamma <- predict(intensity.model, newdata = data[included.obs, ], type = "risk", reference = "zero")
        intensity_weights <- intensity$baseline_intensity * exp_gamma

        # Select term2 computation function once to avoid repeated conditionals
        term2_fn <- if (term2_method == "fast") {
            compute_term2_influence_fast
        } else {
            compute_term2_influence_original
        }

        # Helper function to run optimization for a single alpha value
        fit_single_alpha <- function(current_alpha) {
            # Pre-compute term1 components for this alpha (depends on alpha, not beta)
            expected <- compute_SensIAT_expected_values(
                model = outcome.model,
                alpha = current_alpha,
                new.data = data[included.obs, ]
            )
            term1.deviation.by.observation <-
                (Y[included.obs] - expected$E_Yexp_alphaY / expected$E_exp_alphaY) /
                intensity_weights

            # Build per-patient expected value caches for this alpha
            expected_cache_map <- new.env(parent = emptyenv())
            get_expected_cache_for <- function(patient_id) {
                key <- as.character(patient_id)
                if (!exists(key, expected_cache_map, inherits = FALSE)) {
                    patient_data_local <- data[id == patient_id, ]
                    cache_env <- new.env(parent = emptyenv())
                    
                    # Check if we can use fast pmf-based caching (for single-index models)
                    use_fast_cache <- is(outcome.model, "SensIAT::Single-index-outcome-model") &&
                                      !is.null(attr(outcome.model, "kernel")) &&
                                      !is.null(outcome.model$bandwidth)
                    
                    if (use_fast_cache) {
                        # Pre-compute constants for fast pmf path (same as in make_term2_integrand_fast)
                        Xi <- model.matrix(terms(outcome.model), outcome.model$data)
                        Yi <- model.response(model.frame(outcome.model))
                        beta_outcome <- outcome.model$coef
                        Xb_all <- as.vector(Xi %*% beta_outcome)
                        y_seq <- sort(unique(Yi))
                        kernel <- attr(outcome.model, "kernel")
                        bandwidth <- outcome.model$bandwidth
                        
                        # Pre-extract patient times and outcomes for fast interval lookup
                        outcome_var <- as.character(rlang::f_lhs(formula(outcome.model)))
                        time_candidates <- c("..time..", "Time", "time", "t", "T", "obstime", "obs_time")
                        time_var_local <- NULL
                        for (candidate in time_candidates) {
                            if (candidate %in% names(patient_data_local)) {
                                time_var_local <- candidate
                                break
                            }
                        }
                        if (is.null(time_var_local)) {
                            numeric_cols <- names(patient_data_local)[sapply(patient_data_local, is.numeric)]
                            time_var_local <- setdiff(numeric_cols, outcome_var)[1]
                        }
                        patient_times <- patient_data_local[[time_var_local]]
                        patient_outcomes <- patient_data_local[[outcome_var]]
                        valid_idx <- !is.na(patient_outcomes)
                        patient_times <- patient_times[valid_idx]
                        patient_outcomes <- patient_outcomes[valid_idx]
                        
                        # Prepare model.matrix template for fast xb computation
                        term_spec <- delete.response(terms(outcome.model))
                        template_df <- data.frame(..prev_outcome.. = 0, ..delta_time.. = 0, check.names = FALSE)
                        template_row <- model.matrix(term_spec, data = template_df)
                        mm_colnames <- colnames(template_row)
                        idx_delta <- which(mm_colnames == "..delta_time..")
                        if (length(idx_delta) == 0L) idx_delta <- grep("delta_time", mm_colnames, fixed = TRUE)
                        if (length(idx_delta) != 1L) idx_delta <- NA_integer_
                        
                        # Cache ns-basis rows for unique observed prev_outcome values
                        unique_prev_outcomes <- unique(patient_outcomes)
                        ns_cache <- new.env(parent = emptyenv())
                        for (val in unique_prev_outcomes) {
                            df1 <- data.frame(..prev_outcome.. = val, ..delta_time.. = 0, check.names = FALSE)
                            row0 <- model.matrix(term_spec, data = df1)
                            ns_cache[[as.character(val)]] <- row0
                        }
                        
                        expected_get <- function(t) {
                            k <- as.character(signif(t, 12))
                            if (!exists(k, cache_env, inherits = FALSE)) {
                                # Use fast pmf path (same logic as make_term2_integrand_fast)
                                idx <- findInterval(t, patient_times, left.open = FALSE)
                                if (idx < 1L) idx <- 1L
                                if (idx > length(patient_times)) idx <- length(patient_times)
                                
                                prev_outcome <- patient_outcomes[idx]
                                delta_time <- t - patient_times[idx]
                                
                                if (is.na(prev_outcome) || is.null(prev_outcome)) {
                                    df1 <- data.frame(..prev_outcome.. = 0, ..delta_time.. = delta_time, check.names = FALSE)
                                    x_row <- model.matrix(term_spec, data = df1)[1, , drop = TRUE]
                                } else {
                                    key_po <- as.character(prev_outcome)
                                    if (exists(key_po, envir = ns_cache, inherits = FALSE) && !is.na(idx_delta)) {
                                        x_row <- ns_cache[[key_po]][1, , drop = TRUE]
                                        x_row[idx_delta] <- delta_time
                                    } else {
                                        df1 <- data.frame(..prev_outcome.. = prev_outcome, ..delta_time.. = delta_time, check.names = FALSE)
                                        x_row <- model.matrix(term_spec, data = df1)[1, , drop = TRUE]
                                    }
                                }
                                
                                xb <- sum(x_row * beta_outcome)
                                pmf <- pcoriaccel_estimate_pmf(Xb = Xb_all, Y = Yi, xi = xb, y_seq = y_seq, h = bandwidth, kernel = kernel)
                                
                                E_exp_alphaY <- sum(exp(current_alpha * y_seq) * pmf)
                                E_Yexp_alphaY <- sum(y_seq * exp(current_alpha * y_seq) * pmf)
                                
                                cache_env[[k]] <- list(
                                    E_exp_alphaY = E_exp_alphaY,
                                    E_Yexp_alphaY = E_Yexp_alphaY
                                )
                            }
                            cache_env[[k]]
                        }
                    } else {
                        # Fallback: use generic compute_SensIAT_expected_values for other model types
                        expected_get <- function(t) {
                            k <- as.character(signif(t, 12))
                            if (!exists(k, cache_env, inherits = FALSE)) {
                                df <- impute_data(t, patient_data_local)
                                ev <- compute_SensIAT_expected_values(
                                    model = outcome.model,
                                    alpha = current_alpha,
                                    new.data = df
                                )
                                cache_env[[k]] <- list(
                                    E_exp_alphaY = ev$E_exp_alphaY,
                                    E_Yexp_alphaY = ev$E_Yexp_alphaY
                                )
                            }
                            cache_env[[k]]
                        }
                    }
                    
                    expected_cache_map[[key]] <- expected_get
                }
                expected_cache_map[[key]]
            }

            # Helper function to compute influence for a single patient
            compute_influence_term_1_by_observation <- function(beta) {
                weights <- tryCatch({
                    purrr::map(time[included.obs], W, beta = beta)
                }, error = function(e) {
                    # Return zero weights if weight computation fails
                    replicate(length(time[included.obs]), rep(0, ncol(base)), simplify=FALSE)
                })
                
                # Guard against NaN/NA values
                weights <- purrr::map(weights, function(w) {
                    ifelse(is.finite(w), w, 0)
                })
                
                purrr::map2(weights, term1.deviation.by.observation, function(w, t) {
                    result <- w * t
                    ifelse(is.finite(result), result, 0)
                })
            }
            compute_influence_by_patient <-
            function(patient_id, beta,
                     term1.by.observation = compute_influence_term_1_by_observation(beta)) {
                patient_data <- data[id == patient_id, ]
                expected_get_fn <- if (isTRUE(use_expected_cache)) get_expected_cache_for(patient_id) else NULL
                term2 <- tryCatch({
                    term2_fn(
                        patient_data = patient_data,
                        outcome_model = outcome.model,
                        base = base,
                        alpha = current_alpha,
                        marginal_beta = beta,
                        V_inv = V.inv,
                        tmin = tmin,
                        tmax = tmax,
                        impute_fn = impute_data,
                        inv_link = inv.link,
                        W = W,
                        expected_get = expected_get_fn,
                        time_var = time.var
                    )
                }, error = function(e) {
                    # If term2 computation fails, return zero
                    rep(0, ncol(base))
                })

                # Term1 contributions for this patient (from observations where time[fu] is valid)
                patient_indices_in_fu <- which(id[included.obs] == patient_id)
                if(length(patient_indices_in_fu) == 0){
                    term1 <- rep(0, ncol(base))
                } else {
                    term1 <- reduce(term1.by.observation[patient_indices_in_fu], `+`)
                }

                list(id = patient_id, term1 = term1, term2 = term2, total= term1 + term2)
            }

            # Define the influence function that aggregates by-patient influences
            influence <- function(beta) {
                influence_by_patient <- tryCatch({
                    map(unique_ids, compute_influence_by_patient, beta = beta)
                }, error = function(e) {
                    # If any patient computation fails, return zeros
                    lapply(unique_ids, function(id) list(total = rep(0, ncol(base))))
                })
                
                # Aggregate: sum all term1 observations + sum all term2 patient contributions
                result <- tryCatch({
                    reduce(map(influence_by_patient, getElement, 'total'), `+`)
                }, error = function(e) {
                    rep(0, ncol(base))
                })
                
                # Ensure finite result - element-wise replacement
                result <- purrr::map_dbl(result, function(x) {
                    if (!is.finite(x)) 0 else x
                })
                
                return(result)
            }

            # For identity link, the influence equation is linear in beta and can be solved directly
            # This matches the approach in fit_SensIAT_marginal_mean_model
            if (link == "identity") {
                # Compute influence for all patients (beta doesn't affect W for identity link)
                influence_by_patient <- map(unique_ids, compute_influence_by_patient, beta = rep(0, ncol(base)))
                
                # Extract term1 and term2 separately
                term1_matrix <- do.call(rbind, map(influence_by_patient, getElement, 'term1'))
                term2_matrix <- do.call(rbind, map(influence_by_patient, getElement, 'term2'))
                
                # Closed-form solution: beta_hat = V^{-1} * mean(term1 + term2)
                uncorrected_beta_hat <- (colSums(term1_matrix) + colSums(term2_matrix)) / length(unique_ids)
                V <- GramMatrix(base)
                V_inv <- solve(V)
                solution_par <- as.vector(V_inv %*% uncorrected_beta_hat)
                
                solution <- list(
                    par = solution_par,
                    convergence = 0,
                    message = "Closed-form solution (identity link)"
                )
            } else {
                # Non-identity links require iterative optimization
                solution <- BB::sane(
                    par = rep(0, ncol(base)),
                    fn = influence,
                    control = BBsolve.control
                )
            }
            if (solution$convergence != 0) {
                stop("fit_SensIAT_marginal_mean_model_generalized: Optimization did not converge for alpha=", 
                     current_alpha, ". Message: ", solution$message)
            }

            # Compute final influence by patient at the solution
            final_influence_by_patient <- map(unique_ids, compute_influence_by_patient, beta = solution$par)
            
            # Rebuild term1 and term2 matrices for variance computation
            final_term1_matrix <- do.call(rbind, map(final_influence_by_patient, getElement, 'term1'))
            final_term2_matrix <- do.call(rbind, map(final_influence_by_patient, getElement, 'term2'))
            final_total_matrix <- final_term1_matrix + final_term2_matrix
            
            # Compute coefficient variance using influence-based formula
            # variance = V_inv * cov(influence) * V_inv^T / n
            V <- GramMatrix(base)
            V_inv <- solve(V)
            n_patients <- length(unique_ids)
            
            # Center the influence matrix
            mean_influence <- colMeans(final_total_matrix)
            centered_influence <- sweep(final_total_matrix, 2, mean_influence, "-")
            
            # Compute covariance
            cov_influence <- crossprod(centered_influence) / n_patients
            
            # Variance of coefficients
            coefficient_variance <- V_inv %*% cov_influence %*% t(V_inv) / n_patients
            
            # Create influence tibble
            influence_df <- tibble(
                id = purrr::map(final_influence_by_patient, \(x) x$id) |> purrr::list_c(),
                term1 = purrr::map(final_influence_by_patient, \(x) x$term1),
                term2 = purrr::map(final_influence_by_patient, \(x) x$term2),
                total = purrr::map(final_influence_by_patient, \(x) x$total)
            )
            
            list(
                alpha = current_alpha,
                coefficients = solution$par,
                coefficient_variance = coefficient_variance,
                influence = influence_df
            )
        }

        # Process each alpha value
        results_by_alpha <- purrr::map(alpha, fit_single_alpha)

        # Extract results for all alphas
        coefficients_list <- purrr::map(results_by_alpha, getElement, "coefficients")
        coefficient_variance_list <- purrr::map(results_by_alpha, getElement, "coefficient_variance")
        influence_results <- purrr::map(results_by_alpha, getElement, "influence")
        
        # Always keep as lists for consistency with original function
        coefficients_out <- coefficients_list
        coefficient_variance_out <- coefficient_variance_list
        influence_out <- influence_results

        structure(
            list(
                models = list(
                    intensity = intensity.model,
                    outcome = outcome.model
                ),
                data = data,
                influence = influence_out,
                alpha = alpha,
                coefficients = coefficients_out,
                coefficient.variance = coefficient_variance_out,
                influence.args = list(...),
                base = base,
                V.inverse = V.inv
            ),
            class = "SensIAT_marginal_mean_model_generalized",
            call = match.call(expand.dots = TRUE)
        )
    }
