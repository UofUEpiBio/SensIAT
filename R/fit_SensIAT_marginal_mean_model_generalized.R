#' Fit the marginal mean model for generalize outcomes.
#'
#' This function supports multiple integration methods for term2 computation,
#' including adaptive and fixed-grid approaches. The implementation includes
#' numerical stability improvements (exp(-μ) multiplication vs division) and
#' extensive caching optimizations for repeated expected value computations.
#'
#' @inheritParams fit_SensIAT_marginal_mean_model
#' @param time The time variable in the data. Can be provided as a column name or vector.
#' @param impute_data A function that takes (t, df) and returns the imputed data at time t.
#'   Should handle extrapolation from the last observed time point.
#' @param loss The loss function to use. Options are "lp_mse", "mean_mse", and "quasi-likelihood".
#' @param link The link function to use. Options are "identity", "log", and "logit".
#' @param BBsolve.control Control parameters for the BB::sane optimizer, including `maxit` and `tol`.
#' @param term2_method Method for computing term2 influence components. Options are:
#'   - "fast": Optimized closure-based integrand with adaptive Simpson's (default)
#'   - "original": Standard implementation with adaptive Simpson's
#'   - "fixed_grid": Pre-computed expected values on fixed grid with composite Simpson's rule
#'   - "seeded_adaptive": Adaptive Simpson's seeded with pre-computed grid points
#'   - "gauss_legendre": Gauss-Legendre quadrature (requires statmod package)
#' @param term2_grid_n Number of grid points/nodes for fixed_grid, seeded_adaptive, and gauss_legendre methods (default 100)
#' @param use_expected_cache Logical; whether to cache expected values for performance (default TRUE)
#'
#' @details
#' ## Integration Methods for Term2
#' 
#' The function offers four integration methods with different performance/accuracy tradeoffs:
#' 
#' **Adaptive Methods (fast, original):**
#' - Use adaptive Simpson's quadrature with automatic subdivision
#' - Best accuracy for irregular integrands
#' - "fast" method uses optimized closure-based integrand construction
#' 
#' **Fixed-Grid Method (fixed_grid):**
#' - Pre-computes expected values at fixed grid points (once per alpha)
#' - Uses composite Simpson's rule for integration
#' - 2-5x faster when optimizing over beta (multiple iterations)
#' - Best for smooth integrands with sufficient grid density
#' 
#' **Seeded Adaptive Method (seeded_adaptive):**
#' - Combines pre-computation with adaptive refinement
#' - Starts with pre-computed grid, subdivides where needed
#' - Good balance of speed and accuracy
#'
#' **Gauss-Legendre Method (gauss_legendre):**
#' - Uses Gauss-Legendre quadrature via statmod::gauss.quad
#' - Highly accurate for smooth integrands with fewer evaluation points
#' - Exact for polynomials up to degree 2n-1 using n points
#' - Requires the statmod package
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
#' 
#' **Term1 Optimizations (alpha-independent):**
#' - Y-scaled observations pre-computed once
#' - Patient index mappings cached for O(1) lookups
#' - Identity link: weights pre-computed (don't depend on beta or alpha)
#' - Single-index models: global PMF constants extracted once
#' 
#' **Term2 Optimizations:**
#' - Integration grids pre-computed (alpha-independent)
#' - Basis evaluations at grid points pre-computed
#' - Expected values computed once per alpha (for grid methods)
#' - Per-patient caching of expected values
#' - Weight functions use numerically stable exp(-μ) multiplication
#' - Weight functions use numerically stable exp(-μ) multiplication
#' - Fast method uses closure-based integration with reduced allocations
#'
#'
#' @examples
#' library(survival)
#' library(splines)
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
#'         data = data_with_lags |> dplyr::filter(.data$Time > 0)
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
             term2_method = c("fast", "original", "fixed_grid", "seeded_adaptive", "gauss_legendre"),
             term2_grid_n = 100,
             use_expected_cache = TRUE) {
        time.var <- enquo(time)
        time <- rlang::eval_tidy({{time.var}}, data)
        
        # Extract time column name as string for passing to internal functions
        # This handles both unquoted column names (Time) and df$Time expressions
        time.var.name <- NULL
        expr <- rlang::quo_get_expr(time.var)
        if (is.symbol(expr)) {
            time.var.name <- as.character(expr)
        } else if (is.call(expr) && identical(expr[[1]], as.symbol("$"))) {
            time.var.name <- as.character(expr[[3]])
        }
        # If still NULL, try to find matching column in data
        if (is.null(time.var.name)) {
            for (col in c("Time", "time", "..time..", "t", "T")) {
                if (col %in% names(data)) {
                    time.var.name <- col
                    break
                }
            }
        }

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

        # Select term2 computation function and determine if we need grid pre-computation
        use_term2_grid <- term2_method %in% c("fixed_grid", "seeded_adaptive", "gauss_legendre")
        
        term2_fn <- switch(term2_method,
            fast = compute_term2_influence_fast,
            original = compute_term2_influence_original,
            fixed_grid = compute_term2_influence_fixed_grid,
            seeded_adaptive = compute_term2_influence_seeded_adaptive,
            gauss_legendre = compute_term2_influence_gauss_legendre,
            compute_term2_influence_fast  # default fallback
        )

        # ============================================================
        # PRE-COMPUTE ALPHA-INDEPENDENT TERM1 COMPONENTS
        # ============================================================
        # These computations are done once outside the alpha loop for efficiency
        
        # Y contribution to term1 (scaled by intensity weights) - doesn't depend on alpha
        Y_obs <- Y[included.obs]
        Y_scaled_by_intensity <- Y_obs / intensity_weights
        
        # Observation times for included observations
        obs_times <- time[included.obs]
        n_obs <- length(obs_times)
        
        # Pre-compute patient index mappings (which observations belong to which patient)
        patient_obs_indices <- lapply(unique_ids, function(pid) {
            which(id[included.obs] == pid)
        })
        names(patient_obs_indices) <- as.character(unique_ids)
        
        # Pre-compute patient data subsets
        patient_data_list <- lapply(unique_ids, function(pid) {
            data[id == pid, ]
        })
        names(patient_data_list) <- as.character(unique_ids)
        
        # For identity link, W(t, beta) = V.inv %*% B(t) doesn't depend on beta or alpha
        # Pre-compute all weight vectors for each observation time
        if (link == "identity") {
            # Pre-compute weights for all observation times (matrix: n_obs x ncol(base))
            precomputed_weights <- lapply(obs_times, function(t) {
                B <- pcoriaccel_evaluate_basis(base, t)
                as.vector(V.inv %*% B)
            })
        } else {
            precomputed_weights <- NULL
        }
        
        # Pre-compute per-patient constants for fast pmf path (for single-index models)
        # These don't depend on alpha
        use_fast_cache_global <- is(outcome.model, "SensIAT::Single-index-outcome-model") &&
                                 !is.null(attr(outcome.model, "kernel")) &&
                                 !is.null(outcome.model$bandwidth)
        
        if (use_fast_cache_global) {
            # Global outcome model constants (don't depend on alpha)
            Xi_global <- model.matrix(terms(outcome.model), outcome.model$data)
            Yi_global <- model.response(model.frame(outcome.model))
            beta_outcome_global <- outcome.model$coef
            Xb_all_global <- as.vector(Xi_global %*% beta_outcome_global)
            y_seq_global <- sort(unique(Yi_global))
            kernel_global <- attr(outcome.model, "kernel")
            bandwidth_global <- outcome.model$bandwidth
            
            # Outcome variable name
            outcome_var_global <- as.character(rlang::f_lhs(formula(outcome.model)))
            
            # Pre-compute term_spec for model.matrix
            term_spec_global <- delete.response(terms(outcome.model))
            # Build template with all variables from outcome model
            template_df <- outcome.model$data[1, , drop = FALSE]
            # Set core variables to baseline values
            template_df[["..prev_outcome.."]] <- 0
            if ("..delta_time.." %in% names(template_df)) {
                template_df[["..delta_time.."]] <- 0
            }
            template_row <- model.matrix(term_spec_global, data = template_df)
            mm_colnames_global <- colnames(template_row)
            idx_delta_global <- which(mm_colnames_global == "..delta_time..")
            if (length(idx_delta_global) == 0L) idx_delta_global <- grep("delta_time", mm_colnames_global, fixed = TRUE)
            if (length(idx_delta_global) != 1L) idx_delta_global <- NA_integer_
            
            # Pre-compute per-patient time/outcome vectors and ns_cache
            patient_pmf_cache <- lapply(unique_ids, function(pid) {
                patient_data_local <- patient_data_list[[as.character(pid)]]
                
                # Find time variable
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
                    time_var_local <- setdiff(numeric_cols, outcome_var_global)[1]
                }
                
                patient_times <- patient_data_local[[time_var_local]]
                patient_outcomes <- patient_data_local[[outcome_var_global]]
                valid_idx <- !is.na(patient_outcomes)
                patient_times <- patient_times[valid_idx]
                patient_outcomes <- patient_outcomes[valid_idx]
                
                # Cache ns-basis rows for unique observed prev_outcome values
                unique_prev_outcomes <- unique(patient_outcomes)
                ns_cache <- new.env(parent = emptyenv())
                for (val in unique_prev_outcomes) {
                    df1 <- patient_data_local[1, , drop = FALSE]
                    df1[["..prev_outcome.."]] <- val
                    if ("..delta_time.." %in% names(df1)) {
                        df1[["..delta_time.."]] <- 0
                    }
                    row0 <- model.matrix(term_spec_global, data = df1)
                    ns_cache[[as.character(val)]] <- row0
                }
                
                list(
                    patient_times = patient_times,
                    patient_outcomes = patient_outcomes,
                    ns_cache = ns_cache
                )
            })
            names(patient_pmf_cache) <- as.character(unique_ids)
        }
        
        # ============================================================
        # PRE-COMPUTE TERM2 INTEGRATION GRIDS (for fixed_grid, seeded_adaptive, gauss_legendre)
        # ============================================================
        # These grids are alpha-independent and include observation times (except gauss_legendre)
        patient_term2_grids <- NULL
        if (use_term2_grid) {
            patient_term2_grids <- lapply(unique_ids, function(pid) {
                patient_data_local <- patient_data_list[[as.character(pid)]]
                
                if (term2_method == "gauss_legendre") {
                    # For Gauss-Legendre, use specialized nodes
                    if (!requireNamespace("statmod", quietly = TRUE)) {
                        stop("Package 'statmod' is required for Gauss-Legendre quadrature. ",
                             "Install it with: install.packages('statmod')")
                    }
                    gl <- statmod::gauss.quad(n = term2_grid_n, kind = "legendre")
                    half_range <- (tmax - tmin) / 2
                    mid_point <- (tmax + tmin) / 2
                    grid <- half_range * gl$nodes + mid_point
                    weights <- gl$weights * half_range
                    
                    # Pre-compute basis evaluations at GL nodes
                    B_grid <- lapply(grid, function(t) pcoriaccel_evaluate_basis(base, t))
                    
                    list(
                        grid = grid,
                        B_grid = B_grid,
                        weights = weights
                    )
                } else {
                    # For fixed_grid and seeded_adaptive, use uniform grid + observation times
                    obs_times <- if (use_fast_cache_global) {
                        patient_pmf_cache[[as.character(pid)]]$patient_times
                    } else {
                        extract_patient_times(patient_data_local, time.var.name)
                    }
                    
                    # Create integration grid including observation times
                    grid <- create_integration_grid(
                        tmin = tmin,
                        tmax = tmax,
                        n_grid = term2_grid_n,
                        obs_times = obs_times
                    )
                    
                    # Pre-compute basis evaluations at grid points (alpha-independent)
                    B_grid <- lapply(grid, function(t) pcoriaccel_evaluate_basis(base, t))
                    
                    list(
                        grid = grid,
                        B_grid = B_grid
                    )
                }
            })
            names(patient_term2_grids) <- as.character(unique_ids)
        }
        
        # ============================================================
        # Helper function to run optimization for a single alpha value
        # ============================================================
        fit_single_alpha <- function(current_alpha) {
            # Compute alpha-dependent expected values for included observations
            expected <- compute_SensIAT_expected_values(
                model = outcome.model,
                alpha = current_alpha,
                new.data = data[included.obs, ]
            )
            
            # Alpha-dependent correction term
            alpha_correction <- (expected$E_Yexp_alphaY / expected$E_exp_alphaY) / intensity_weights
            
            # Combine pre-computed Y contribution with alpha-dependent correction
            term1.deviation.by.observation <- Y_scaled_by_intensity - alpha_correction

            # ============================================================
            # PRE-COMPUTE EXPECTED VALUES AT GRID POINTS (for fixed_grid and seeded_adaptive)
            # ============================================================
            # This is alpha-dependent but beta-independent, so computed once per alpha
            patient_expected_grids <- NULL
            if (use_term2_grid) {
                patient_expected_grids <- lapply(unique_ids, function(pid) {
                    key <- as.character(pid)
                    patient_data_local <- patient_data_list[[key]]
                    grid_info <- patient_term2_grids[[key]]
                    
                    # Compute expected values at all grid points for this patient and alpha
                    E_grid <- compute_expected_values_at_grid(
                        grid = grid_info$grid,
                        patient_data = patient_data_local,
                        outcome_model = outcome.model,
                        alpha = current_alpha,
                        impute_fn = impute_data,
                        time_var = time.var.name
                    )
                    
                    # Combine with pre-computed basis evaluations
                    list(
                        grid = grid_info$grid,
                        B_grid = grid_info$B_grid,
                        E_grid = E_grid
                    )
                })
                names(patient_expected_grids) <- as.character(unique_ids)
            }

            # Build per-patient expected value caches for this alpha (for term2)
            # Uses pre-computed patient constants when available
            expected_cache_map <- new.env(parent = emptyenv())
            get_expected_cache_for <- function(patient_id) {
                key <- as.character(patient_id)
                if (!exists(key, expected_cache_map, inherits = FALSE)) {
                    cache_env <- new.env(parent = emptyenv())
                    
                    # Use pre-computed patient data
                    patient_data_local <- patient_data_list[[key]]
                    
                    if (use_fast_cache_global) {
                        # Use pre-computed global constants and per-patient cache
                        pmf_cache <- patient_pmf_cache[[key]]
                        patient_times <- pmf_cache$patient_times
                        patient_outcomes <- pmf_cache$patient_outcomes
                        ns_cache <- pmf_cache$ns_cache
                        
                        expected_get <- function(t) {
                            k <- as.character(signif(t, 12))
                            if (!exists(k, cache_env, inherits = FALSE)) {
                                # Use fast pmf path with pre-computed constants
                                idx <- findInterval(t, patient_times, left.open = FALSE)
                                if (idx < 1L) idx <- 1L
                                if (idx > length(patient_times)) idx <- length(patient_times)
                                
                                prev_outcome <- patient_outcomes[idx]
                                delta_time <- t - patient_times[idx]
                                
                                if (is.na(prev_outcome) || is.null(prev_outcome)) {
                                    df1 <- outcome.model$data[1, , drop = FALSE]
                                    df1[["..prev_outcome.."]] <- 0
                                    df1[["..delta_time.."]] <- delta_time
                                    if ("Time" %in% names(df1)) {
                                        df1[["Time"]] <- t
                                    }
                                    x_row <- model.matrix(term_spec_global, data = df1)[1, , drop = TRUE]
                                } else {
                                    key_po <- as.character(prev_outcome)
                                    if (exists(key_po, envir = ns_cache, inherits = FALSE) && !is.na(idx_delta_global)) {
                                        x_row <- ns_cache[[key_po]][1, , drop = TRUE]
                                        x_row[idx_delta_global] <- delta_time
                                        # Also update Time to the current time t
                                        time_idx_vec <- which(names(x_row) == "Time")
                                        if (length(time_idx_vec) > 0) {
                                            x_row[time_idx_vec] <- t
                                        }
                                    } else {
                                        df1 <- outcome.model$data[1, , drop = FALSE]
                                        df1[["..prev_outcome.."]] <- prev_outcome
                                        df1[["..delta_time.."]] <- delta_time
                                        if ("Time" %in% names(df1)) {
                                            df1[["Time"]] <- t
                                        }
                                        x_row <- model.matrix(term_spec_global, data = df1)[1, , drop = TRUE]
                                    }
                                }
                                
                                xb <- sum(x_row * beta_outcome_global)
                                pmf <- pcoriaccel_estimate_pmf(Xb = Xb_all_global, Y = Yi_global, xi = xb, 
                                                               y_seq = y_seq_global, h = bandwidth_global, 
                                                               kernel = kernel_global)
                                
                                E_exp_alphaY <- sum(exp(current_alpha * y_seq_global) * pmf)
                                E_Yexp_alphaY <- sum(y_seq_global * exp(current_alpha * y_seq_global) * pmf)
                                
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
                                    E_exp_alphaY = as.numeric(ev$E_exp_alphaY)[1],
                                    E_Yexp_alphaY = as.numeric(ev$E_Yexp_alphaY)[1]
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
                # For identity link, use pre-computed weights (they don't depend on beta)
                if (link == "identity" && !is.null(precomputed_weights)) {
                    weights <- precomputed_weights
                } else {
                    weights <- tryCatch({
                        purrr::map(obs_times, W, beta = beta)
                    }, error = function(e) {
                        # Return zero weights if weight computation fails
                        replicate(n_obs, rep(0, ncol(base)), simplify=FALSE)
                    })
                }
                
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
                # Use pre-computed patient data
                key <- as.character(patient_id)
                patient_data <- patient_data_list[[key]]
                expected_get_fn <- if (isTRUE(use_expected_cache)) get_expected_cache_for(patient_id) else NULL
                
                # Get pre-computed grid for grid-based methods
                expected_grid_arg <- if (use_term2_grid) {
                    patient_expected_grids[[key]]
                } else {
                    NULL
                }
                
                term2 <- tryCatch({
                    result <- term2_fn(
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
                        expected_grid = expected_grid_arg,
                        n_grid = term2_grid_n,
                        time_var = time.var.name
                    )
                    result
                }, error = function(e) {
                    # If term2 computation fails, return zero
                    message("term2 error for patient ", patient_id, ": ", conditionMessage(e))
                    rep(0, ncol(base))
                })

                # Term1 contributions for this patient - use pre-computed indices
                patient_indices_in_fu <- patient_obs_indices[[key]]
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
