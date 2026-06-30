## Parametric bootstrap orchestration and helpers

#' Sample parametric coefficients from a fitted model using asymptotic multivariate normal distribution.
#'
#' @param model A fitted model with `coef()` and `vcov()` methods
#' @return Named numeric vector of sampled coefficients
sample_parametric_coeffs <- function(model) {
    coefs <- tryCatch(coef(model), error = function(e) NULL)
    if (is.null(coefs)) stop("Model does not provide coef()")
    vc <- tryCatch(vcov(model), error = function(e) NULL)
    if (is.null(vc)) stop("Model does not provide vcov(); cannot sample coefficients")

    if (inherits(model, "SensIAT::Single-index-outcome-model") && is.matrix(vc) && nrow(vc) == length(coefs)) {
        # If the first coefficient is fixed by normalization, sample the free coefficients only.
        if (isTRUE(all.equal(coefs[1], 1)) && all(vc[1, -1] == 0) && all(vc[-1, 1] == 0)) {
            samp_tail <- as.numeric(MASS::mvrnorm(1, mu = as.numeric(coefs[-1]), Sigma = vc[-1, -1]))
            samp <- c(coefs[1], samp_tail)
            names(samp) <- names(coefs)
            return(samp)
        }
    }

    # Use MASS::mvrnorm
    samp <- as.numeric(MASS::mvrnorm(1, mu = as.numeric(coefs), Sigma = vc))
    names(samp) <- names(coefs)
    return(samp)
}


#' Get coefficients for bootstrap simulation.
#'
#' By default this returns original fitted coefficients. If sampling is requested,
#' coefficients are sampled from an asymptotic multivariate normal distribution only when `vcov()` is available.
#'
#' @param model A fitted model with `coef()` and `vcov()` methods.
#' @param sample_coefficients Logical; if `TRUE`, sample coefficients from an asymptotic multivariate normal distribution when `vcov()` is available. If `FALSE` (default), use original fitted coefficients.
#' @param model_label Character; label for the model used in warning messages.
get_bootstrap_coeffs <- function(model, sample_coefficients = FALSE, model_label = "model") {
    coefs <- tryCatch(coef(model), error = function(e) NULL)
    if (is.null(coefs)) stop("Model does not provide coef()")

    if (!isTRUE(sample_coefficients)) {
        return(coefs)
    }

    vc <- tryCatch(vcov(model), error = function(e) NULL)
    if (is.null(vc)) {
        warning(
            sprintf(
                "Sampling requested for %s, but vcov() is unavailable; using original coefficients.",
                model_label
            )
        )
        return(coefs)
    }

    tryCatch(
        sample_parametric_coeffs(model),
        error = function(e) {
            warning(
                sprintf(
                    "Could not sample coefficients for %s (%s); using original coefficients.",
                    model_label,
                    conditionMessage(e)
                )
            )
            coefs
        }
    )
}


normalize_bootstrap_verbosity <- function(verbosity = c("none", "basic", "detailed"), verbose = NULL) {
    if (!is.null(verbose) && isTRUE(verbose)) {
        return("detailed")
    }

    if (is.logical(verbosity) && length(verbosity) == 1) {
        return(if (isTRUE(verbosity)) "detailed" else "none")
    }

    if (is.numeric(verbosity) && length(verbosity) == 1) {
        if (verbosity <= 0) return("none")
        if (verbosity == 1) return("basic")
        return("detailed")
    }

    match.arg(verbosity)
}


bootstrap_log <- function(verbosity, level = c("basic", "detailed"), ...) {
    level <- match.arg(level)
    current_level <- switch(verbosity, none = 0L, basic = 1L, detailed = 2L, 0L)
    needed_level <- if (identical(level, "basic")) 1L else 2L
    if (current_level >= needed_level) {
        cat("[SensIAT bootstrap] ", paste0(..., collapse = ""), "\n", sep = "")
    }
}


#' Build a parametric intensity simulator using sampled coefficients
#'
#' Supports `coxph` objects: sampled coefficients are used to compute linear predictor
#' while baseline hazard is taken from `survival::basehaz()`.
#'
#' @param intensity_model A fitted model (e.g., `coxph`) or a function(t, prev_outcome, visit_num).
#' @param sampled_coef Named numeric vector of sampled coefficients (from `get_bootstrap_coeffs()`).
#' @param covariate_mapping Optional named character vector mapping expected covariate names to model variable names.
#' @return A function of signature `function(t, prev_outcome, visit_num)` returning non-negative numeric intensity.
#' 
#' `r lifecycle::badge("experimental")`
make_parametric_intensity_simulator <- function(intensity_model, sampled_coef, covariate_mapping = NULL) {
    if (!inherits(intensity_model, "coxph")) stop("Only coxph intensity_model supported for parametric sampling currently")

    model_terms <- delete.response(terms(intensity_model))
    # survival::basehaz() can emit a harmless warning on very thin strata
    # (e.g., single event with min(diff(time)) on empty vector).
    bh <- tryCatch(suppressWarnings(survival::basehaz(intensity_model, centered = FALSE)), error = function(e) NULL)
    if (is.null(bh) || nrow(bh) == 0) {
        warning("Could not extract baseline hazard from coxph; using constant hazard 0.01")
        return(function(t, prev_outcome, visit_num) 0.01)
    }
    
    # Extract and clean times/hazard vectors
    bh_times <- bh$time[!is.na(bh$time)]
    bh_cumhaz <- if (!is.null(bh$hazard)) bh$hazard[!is.na(bh$time)] else bh$haz[!is.na(bh$time)]
    
    # Sort by times
    sort_idx <- order(bh_times)
    bh_times <- bh_times[sort_idx]
    bh_cumhaz <- bh_cumhaz[sort_idx]
    
    times <- c(0, bh_times)
    cumhaz <- c(0, bh_cumhaz)
    if (length(times) < 2) {
        warning("Baseline hazard from coxph is degenerate; using constant hazard 0.01")
        return(function(t, prev_outcome, visit_num) 0.01)
    }
    interval_lengths <- diff(times)
    base_hazard <- diff(cumhaz) / interval_lengths
    base_hazard[base_hazard < 0] <- 0

    baseline_hazard_fn <- function(t) {
        if (!is.numeric(t) || length(t) != 1 || is.na(t)) {
            stop("Time must be a single numeric value")
        }
        if (t < 0) {
            return(0)
        }
        idx <- findInterval(t, times, rightmost.closed = FALSE)
        if (idx == 0 || idx > length(base_hazard)) {
            return(0)
        }
        return(base_hazard[idx])
    }

    fn <- function(t, prev_outcome, visit_num) {
        nd <- data.frame(prev_outcome = prev_outcome, visit_num = visit_num)
        # Build model matrix and compute lp using sampled_coef
        X_new <- tryCatch(model.matrix(model_terms, data = nd), error = function(e) NULL)
        if (is.null(X_new)) {
            lp <- 0
        } else {
            # Align coefficients
            cols <- colnames(X_new)
            coef_vec <- sampled_coef
            missing <- setdiff(cols, names(coef_vec))
            if (length(missing) > 0) {
                coef_vec[missing] <- 0
            }
            lp <- as.numeric(X_new %*% coef_vec[cols])
        }
        h <- baseline_hazard_fn(t)
        val <- h * exp(lp)
        if (!is.finite(val) || val < 0) val <- 0
        return(as.numeric(val))
    }
    attr(fn, "baseline_event_times") <- times[-1]
    attr(fn, "baseline_hazard_segments") <- base_hazard
    return(fn)
}


#' Build a parametric outcome simulator from a fitted single-index outcome model
#' using sampled coefficients for the single-index projection.
#' 
#' @param outcome_model A fitted `SensIAT::Single-index-outcome-model` object.
#' @param sampled_coef Named numeric vector of sampled coefficients (from `get_bootstrap_coeffs()`).
#' @param covariate_mapping Optional named character vector mapping expected covariate names to model variable names.
#' `r lifecycle::badge("experimental")` 
make_parametric_single_index_simulator <- function(outcome_model, sampled_coef, covariate_mapping = NULL) {
    if (!inherits(outcome_model, "SensIAT::Single-index-outcome-model")) stop("outcome_model must be a fitted single-index outcome model")
    bandwidth <- if (!is.null(outcome_model$bandwidth)) outcome_model$bandwidth else stop("bandwidth missing in outcome_model")
    data_orig <- outcome_model$data
    default_names <- c(prev_outcome = "..prev_outcome..", time = "..time..", delta_time = "..delta_time..")
    if (!is.null(covariate_mapping)) {
        if (!is.character(covariate_mapping) || is.null(names(covariate_mapping))) stop("covariate_mapping must be a named character vector")
        default_names[names(covariate_mapping)] <- covariate_mapping
    }

    model_terms <- delete.response(terms(outcome_model))
    y_resp <- model.response(model.frame(outcome_model))
    y_seq <- sort(unique(y_resp))
    kernel_type <- attr(outcome_model, "kernel")
    
    # OPTIMIZATION 1: Pre-compute X_orig and lp0 once (outside the closure)
    # These are constant across all simulator calls since they only depend on:
    # - Original training data (fixed)
    # - Sampled coefficients (fixed for this simulator instance)
    X_orig <- model.matrix(model_terms, data = data_orig)
    lp0 <- as.vector(X_orig %*% as.numeric(sampled_coef[colnames(X_orig)]))

    function(prev_outcome, time, delta_time, newdata = NULL) {
        if (!is.null(newdata)) {
            if (!is.data.frame(newdata) || nrow(newdata) < 1) stop("newdata must be a data.frame with at least one row")
            nd <- newdata[1, , drop = FALSE]
        } else {
            nd <- data.frame(prev_outcome = prev_outcome, time = time, delta_time = delta_time, stringsAsFactors = FALSE)
            names(nd) <- unname(default_names)
        }

        X_new <- tryCatch(model.matrix(model_terms, data = nd), error = function(e) stop("Failed to build model matrix for outcome newdata: ", conditionMessage(e)))
        xb <- as.vector(X_new %*% as.numeric(sampled_coef[colnames(X_new)]))

        # lp0 is pre-computed; only X_new and xb are computed per call
        Fhat <- pcoriaccel_NW(Xb = lp0, Y = y_resp, xb = xb[1], y_seq = y_seq, h = bandwidth, kernel = kernel_type)
        Fhat_vec <- as.vector(Fhat)
        pmf <- pmax(0, diff(c(0, Fhat_vec)))
        if (sum(pmf) <= 0) return(mean(y_resp, na.rm = TRUE))
        pmf <- pmf / sum(pmf)
        sampled <- sample(y_seq, size = 1, prob = pmf)
        return(sampled)
    }
}


#' Parametric bootstrap orchestration
#' `r lifecycle::badge("experimental")`
#' 
#' @param nboot Number of bootstrap replicates
#' @param intensity_model Fitted intensity model ([coxph]) or function
#' @param outcome_model Fitted outcome model (single-index) or NULL
#' @param simulate_args List of arguments to pass to `simulate_SensIAT_data()` (e.g., n_subjects, End, intensity_bound)
#' @param seed Optional seed for reproducibility
#' @param progress Logical; show progress bar when available.
#' @param sample_coefficients Logical; if `TRUE`, sample coefficients from an
#'        asymptotic multivariate normal distribution when `vcov()` is available. If `FALSE` (default), use
#'        original fitted coefficients.
#' @param verbosity Logging verbosity for bootstrap internals: one of
#'        `"none"`, `"basic"`, or `"detailed"`.
#' @param verbose Deprecated shortcut; if `TRUE`, equivalent to
#'        `verbosity = "detailed"`.
#' @return A list of simulated datasets (length `nboot`)
parametric_bootstrap <- function(nboot = 100,
                                 intensity_model = NULL,
                                 outcome_model = NULL,
                                 simulate_args = list(),
                                 seed = NULL,
                                 progress = interactive(),
                                 sample_coefficients = FALSE,
                                 verbosity = c("none", "basic", "detailed"),
                                 verbose = NULL) {
    verbosity <- normalize_bootstrap_verbosity(verbosity = verbosity, verbose = verbose)
    bootstrap_log(
        verbosity,
        "basic",
        "parametric_bootstrap started: nboot=", nboot,
        ", sample_coefficients=", sample_coefficients
    )

    if (!is.null(seed)) set.seed(seed)
    if (!is.null(seed)) {
        bootstrap_log(verbosity, "basic", "seed set to ", seed)
    }
    results <- vector("list", nboot)
    intensity_sim_fixed <- NULL
    outcome_sim_fixed <- NULL
    intensity_bound_fixed <- NULL
    
    if (progress && rlang::is_installed("progress")) {
        pb <- progress::progress_bar$new(
            format = "  parametric bootstrap [:bar] :current/:total(:percent) eta: :eta",
            total = nboot
        )
        pb$tick(0)
        on.exit(pb$terminate(), add = TRUE)
        bootstrap_log(verbosity, "basic", "progress bar enabled for simulation loop")
    }

    # When coefficients are fixed, build simulators once and reuse across replicates.
    if (!sample_coefficients) {
        if (!is.null(intensity_model)) {
            bootstrap_log(verbosity, "detailed", "building fixed intensity simulator from original coefficients")
            sc_fixed <- get_bootstrap_coeffs(
                intensity_model,
                sample_coefficients = FALSE,
                model_label = "intensity model"
            )
            intensity_sim_fixed <- make_parametric_intensity_simulator(intensity_model, sc_fixed)
        }
        if (!is.null(outcome_model)) {
            bootstrap_log(verbosity, "detailed", "building fixed outcome simulator from original coefficients")
            sc2_fixed <- get_bootstrap_coeffs(
                outcome_model,
                sample_coefficients = FALSE,
                model_label = "outcome model"
            )
            outcome_sim_fixed <- make_parametric_single_index_simulator(outcome_model, sc2_fixed)
        }

        if (!is.null(intensity_sim_fixed) && is.null(simulate_args$intensity_bound)) {
            bootstrap_log(verbosity, "detailed", "computing shared intensity_bound from fixed intensity simulator")
            End_tmp <- if (!is.null(simulate_args$End)) simulate_args$End else 1
            max_vis_tmp <- if (!is.null(simulate_args$max_visits)) simulate_args$max_visits else 5
            init_mean_tmp <- if (!is.null(simulate_args$initial_outcome_mean)) simulate_args$initial_outcome_mean else 0
            outcome_sd_tmp <- if (!is.null(simulate_args$initial_outcome_sd)) simulate_args$initial_outcome_sd else 1

            times_grid <- unique(c(
                seq(0, End_tmp, length.out = max(500, ceiling(End_tmp) * 100)),
                sort(c(0, End_tmp, seq(0, End_tmp, length.out = max(50, ceiling(End_tmp) * 10))))
            ))
            prev_outcomes <- unique(c(
                init_mean_tmp,
                0,
                init_mean_tmp + 3 * outcome_sd_tmp,
                init_mean_tmp - 3 * outcome_sd_tmp
            ))
            prev_outcomes <- prev_outcomes[is.finite(prev_outcomes)]
            evals <- double(0)
            baseline_times <- attr(intensity_sim_fixed, "baseline_event_times")
            if (!is.null(baseline_times)) {
                times_grid <- unique(c(times_grid, baseline_times, pmax(0, baseline_times - 1e-4), pmin(End_tmp, baseline_times + 1e-4)))
            }
            for (tt in times_grid) {
                for (pv in prev_outcomes) {
                    for (vn in seq_len(min(5, max_vis_tmp))) {
                        val <- tryCatch(intensity_sim_fixed(tt, pv, vn), error = function(e) 0)
                        if (!is.finite(val) || val < 0) val <- 0
                        evals <- c(evals, val)
                    }
                }
            }
            intensity_bound_fixed <- max(1e-6, max(evals, na.rm = TRUE) * 10000, 1e5)
            bootstrap_log(verbosity, "detailed", "shared intensity_bound computed: ", signif(intensity_bound_fixed, 6))
        }
    }

    for (b in seq_len(nboot)) {
        t_rep <- proc.time()
        bootstrap_log(verbosity, "detailed", "replicate ", b, "/", nboot, ": simulation start")

        # Sample coefficients and build simulators as available
        intensity_sim <- intensity_sim_fixed
        outcome_sim <- outcome_sim_fixed

        if (is.null(intensity_sim) && !is.null(intensity_model)) {
            bootstrap_log(verbosity, "detailed", "replicate ", b, ": sampling intensity coefficients")
            sc <- get_bootstrap_coeffs(
                intensity_model,
                sample_coefficients = sample_coefficients,
                model_label = "intensity model"
            )
            intensity_sim <- make_parametric_intensity_simulator(intensity_model, sc)
        }

        if (is.null(outcome_sim) && !is.null(outcome_model)) {
            bootstrap_log(verbosity, "detailed", "replicate ", b, ": sampling outcome coefficients")
            sc2 <- get_bootstrap_coeffs(
                outcome_model,
                sample_coefficients = sample_coefficients,
                model_label = "outcome model"
            )
            outcome_sim <- make_parametric_single_index_simulator(outcome_model, sc2)
        }

        # Merge supplied simulate_args and inject simulators
        args <- simulate_args
        if (!is.null(intensity_sim)) {
            args$intensity_fn <- intensity_sim
            if (!is.null(intensity_bound_fixed)) {
                args$intensity_bound <- intensity_bound_fixed
            }
            # require an intensity_bound; try to set from simulate_args or compute a safe bound
            if (is.null(args$intensity_bound)) {
                bootstrap_log(verbosity, "detailed", "replicate ", b, ": computing intensity_bound")
                End_tmp <- if (!is.null(simulate_args$End)) simulate_args$End else 1
                max_vis_tmp <- if (!is.null(simulate_args$max_visits)) simulate_args$max_visits else 5
                init_mean_tmp <- if (!is.null(simulate_args$initial_outcome_mean)) simulate_args$initial_outcome_mean else 0
                outcome_sd_tmp <- if (!is.null(simulate_args$initial_outcome_sd)) simulate_args$initial_outcome_sd else 1

                times_grid <- unique(c(
                    seq(0, End_tmp, length.out = max(500, ceiling(End_tmp) * 100)),
                    sort(c(0, End_tmp, seq(0, End_tmp, length.out = max(50, ceiling(End_tmp) * 10))))
                ))
                prev_outcomes <- unique(c(
                    init_mean_tmp,
                    0,
                    init_mean_tmp + 3 * outcome_sd_tmp,
                    init_mean_tmp - 3 * outcome_sd_tmp
                ))
                prev_outcomes <- prev_outcomes[is.finite(prev_outcomes)]
                evals <- double(0)
                baseline_times <- attr(intensity_sim, "baseline_event_times")
                if (!is.null(baseline_times)) {
                    times_grid <- unique(c(times_grid, baseline_times, pmax(0, baseline_times - 1e-4), pmin(End_tmp, baseline_times + 1e-4)))
                }
                for (tt in times_grid) {
                    for (pv in prev_outcomes) {
                        for (vn in seq_len(min(5, max_vis_tmp))) {
                            val <- tryCatch(intensity_sim(tt, pv, vn), error = function(e) 0)
                            if (!is.finite(val) || val < 0) val <- 0
                            evals <- c(evals, val)
                        }
                    }
                }
                upper <- max(1e-6, max(evals, na.rm = TRUE) * 10000, 1e5)
                args$intensity_bound <- upper
                bootstrap_log(verbosity, "detailed", "replicate ", b, ": intensity_bound=", signif(upper, 6))
            }
        }
        if (!is.null(outcome_sim)) args$outcome_simulator <- outcome_sim

        # Call simulation
        sim <- do.call(simulate_SensIAT_data, args)
        results[[b]] <- sim

        t_rep_elapsed <- (proc.time() - t_rep)[[3]]
        bootstrap_log(
            verbosity,
            "detailed",
            "replicate ", b, "/", nboot,
            ": simulation done (rows=", nrow(sim), ", elapsed=", sprintf("%.2f", t_rep_elapsed), " s)"
        )
        
        if (progress && exists("pb")) {
            try(pb$tick(), silent = TRUE)
        }
    }
    bootstrap_log(verbosity, "basic", "parametric_bootstrap completed: ", nboot, " datasets generated")
    return(results)
}

#' Parametric bootstrap for a within-group SensIAT model
#' `r lifecycle::badge("experimental")`
#'
#' @param nboot Number of bootstrap replicates.
#' @param within_group_model A fitted `SensIAT_within_group_model` object.
#' @param simulate_args List of arguments to pass to `simulate_SensIAT_data()`.
#'        If not specified, `End`, `n_subjects`, `initial_outcome_mean`, and
#'        `initial_outcome_sd` are inferred from the fitted model.
#' @param seed Optional seed for reproducibility.
#' @param progress Logical; show progress bar when available.
#' @param sample_coefficients Logical; if `TRUE`, sample coefficients from an
#'        asymptotic multivariate normal distribution when `vcov()` is available. If `FALSE` (default), use
#'        original fitted coefficients.
#' @param refit Logical; if `TRUE` (default), fit a `SensIAT_within_group_model`
#'        on each simulated replicate using the original model's settings.
#' @param return One of `"coefficients"` (default), `"models"`, or `"data"`.
#'        `"coefficients"` is memory-efficient and stores only replicated
#'        marginal mean coefficients.
#' @param prune_models Logical; when `return = "models"`, prune each replicated
#'        model before returning.
#' @param gc_every Integer. Run `gc(FALSE)` every `gc_every` replications.
#'        Use `NULL` to disable explicit garbage collection.
#' @param verbosity Logging verbosity for bootstrap internals: one of
#'        `"none"`, `"basic"`, or `"detailed"`.
#' @param verbose Deprecated shortcut; if `TRUE`, equivalent to
#'        `verbosity = "detailed"`.
#' @return If `return = "coefficients"`, a
#'         `SensIAT_withingroup_bootstrap_results` object. Otherwise returns a
#'         list of replicated fitted models (`"models"`) or simulated datasets
#'         (`"data"`).
#' @examples
#' \dontrun{
#' data("SensIAT_example_data", package = "SensIAT")
#'
#' # Fit a single-index outcome model on a small subset of the example data.
#' small_data <- dplyr::filter(
#'   SensIAT_example_data,
#'   Subject_ID %in% head(unique(SensIAT_example_data$Subject_ID), 8)
#' )
#'
#' model <- fit_SensIAT_within_group_model(
#'   group.data = small_data,
#'   outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
#'   alpha = 0,
#'   id = Subject_ID,
#'   outcome = Outcome,
#'   time = Time,
#'   End = 830,
#'   knots = c(60, 260, 460)
#' )
#'
#' # This example may take a long time because it fits a single-index outcome model
#' # and generates bootstrap replicates.
#' res <- parametric_bootstrap_within_group(
#'   nboot = 2,
#'   within_group_model = model,
#'   simulate_args = list(
#'     n_subjects = 3,
#'     End = 5,
#'     max_visits = 5
#'   ),
#'   seed = 123,
#'   refit = TRUE
#' )
#'
#' print(res[[1]])
#' }
#' @export
parametric_bootstrap_within_group <- function(within_group_model,
                                              nboot = 100,
                                              simulate_args = list(),
                                              seed = NULL,
                                              progress = interactive(),
                                              sample_coefficients = FALSE,
                                              refit = TRUE,
                                              return = c("coefficients", "models", "data"),
                                              prune_models = FALSE,
                                              gc_every = 10L,
                                              verbosity = c("none", "basic", "detailed"),
                                              verbose = NULL) {
    verbosity <- normalize_bootstrap_verbosity(verbosity = verbosity, verbose = verbose)
    return <- match.arg(return)

    bootstrap_log(
        verbosity,
        "basic",
        "parametric_bootstrap_within_group started: nboot=", nboot,
        ", return=", return,
        ", refit=", refit,
        ", sample_coefficients=", sample_coefficients
    )

    if (!inherits(within_group_model, "SensIAT_within_group_model")) {
        stop("within_group_model must be a SensIAT_within_group_model object")
    }

    intensity_model <- within_group_model$models$intensity
    outcome_model <- within_group_model$models$outcome
    model_data <- within_group_model$data

    if (is.null(simulate_args$End)) {
        simulate_args$End <- within_group_model$End
    }
    if (is.null(simulate_args$n_subjects)) {
        if (!is.null(model_data$..id..)) {
            simulate_args$n_subjects <- length(unique(model_data$..id..))
        } else {
            simulate_args$n_subjects <- NA_integer_
        }
    }
    if (is.null(simulate_args$initial_outcome_mean)) {
        baseline_outcomes <- model_data$..outcome..[model_data$..time.. == 0]
        if (length(baseline_outcomes) > 0) {
            simulate_args$initial_outcome_mean <- mean(baseline_outcomes, na.rm = TRUE)
        } else {
            simulate_args$initial_outcome_mean <- mean(model_data$..outcome.., na.rm = TRUE)
        }
    }
    if (is.null(simulate_args$initial_outcome_sd)) {
        baseline_outcomes <- model_data$..outcome..[model_data$..time.. == 0]
        if (length(baseline_outcomes) > 1) {
            simulate_args$initial_outcome_sd <- stats::sd(baseline_outcomes, na.rm = TRUE)
        } else {
            simulate_args$initial_outcome_sd <- stats::sd(model_data$..outcome.., na.rm = TRUE)
        }
    }

    if (!isTRUE(refit) && return != "data") {
        stop("When refit = FALSE, return must be 'data'.")
    }

    if (return == "data") {
        bootstrap_log(verbosity, "basic", "data-only mode selected; delegating to parametric_bootstrap")
        return(parametric_bootstrap(
            nboot = nboot,
            intensity_model = intensity_model,
            outcome_model = outcome_model,
            simulate_args = simulate_args,
            seed = seed,
            progress = progress,
            sample_coefficients = sample_coefficients,
            verbosity = verbosity
        ))
    }

    if (!is.null(seed)) set.seed(seed)
    if (!is.null(seed)) {
        bootstrap_log(verbosity, "basic", "seed set to ", seed)
    }

    intensity_sim_fixed <- NULL
    outcome_sim_fixed <- NULL
    intensity_bound_fixed <- NULL

    compute_safe_intensity_bound <- function(intensity_sim, simulate_args) {
        bootstrap_log(verbosity, "detailed", "computing safe intensity bound")
        End_tmp <- if (!is.null(simulate_args$End)) simulate_args$End else 1
        max_vis_tmp <- if (!is.null(simulate_args$max_visits)) simulate_args$max_visits else 5
        init_mean_tmp <- if (!is.null(simulate_args$initial_outcome_mean)) simulate_args$initial_outcome_mean else 0
        outcome_sd_tmp <- if (!is.null(simulate_args$initial_outcome_sd)) simulate_args$initial_outcome_sd else 1

        # COARSE GRID: 50 time points instead of 500+
        times_grid <- unique(c(
            seq(0, End_tmp, length.out = 50),
            0, End_tmp
        ))
        
        # 5 evenly spaced prev_outcomes
        prev_outcomes <- seq(
            init_mean_tmp - 3 * outcome_sd_tmp,
            init_mean_tmp + 3 * outcome_sd_tmp,
            length.out = 5
        )
        prev_outcomes <- prev_outcomes[is.finite(prev_outcomes)]
        
        # Pre-allocate vector instead of growing it
        evals <- numeric(length(times_grid) * length(prev_outcomes) * 5)
        idx <- 0L
        
        baseline_times <- attr(intensity_sim, "baseline_event_times")
        if (!is.null(baseline_times)) {
            times_grid <- unique(c(times_grid, baseline_times))
        }
        
        # Triple nested loop with pre-allocated storage
        for (tt in times_grid) {
            for (pv in prev_outcomes) {
                for (vn in seq_len(min(5, max_vis_tmp))) {
                    idx <- idx + 1L
                    val <- tryCatch(intensity_sim(tt, pv, vn), error = function(e) 0)
                    if (!is.finite(val) || val < 0) val <- 0
                    evals[idx] <- val
                }
            }
        }
        evals <- evals[seq_len(idx)]
        
        bound <- max(1e-6, max(evals, na.rm = TRUE) * 10000, 1e5)
        bootstrap_log(verbosity, "detailed", "safe intensity bound computed: ", signif(bound, 6))
        bound
    }

    if (!sample_coefficients) {
        if (!is.null(intensity_model)) {
            bootstrap_log(verbosity, "detailed", "building fixed intensity simulator for refit bootstrap")
            sc_fixed <- get_bootstrap_coeffs(
                intensity_model,
                sample_coefficients = FALSE,
                model_label = "intensity model"
            )
            intensity_sim_fixed <- make_parametric_intensity_simulator(intensity_model, sc_fixed)
        }
        if (!is.null(outcome_model)) {
            bootstrap_log(verbosity, "detailed", "building fixed outcome simulator for refit bootstrap")
            sc2_fixed <- get_bootstrap_coeffs(
                outcome_model,
                sample_coefficients = FALSE,
                model_label = "outcome model"
            )
            outcome_sim_fixed <- make_parametric_single_index_simulator(outcome_model, sc2_fixed)
        }
        if (!is.null(intensity_sim_fixed) && is.null(simulate_args$intensity_bound)) {
            intensity_bound_fixed <- compute_safe_intensity_bound(intensity_sim_fixed, simulate_args)
        }
    }

    if (progress && rlang::is_installed("progress")) {
        pb_refit <- progress::progress_bar$new(
            format = "  bootstrap refit [:bar] :current/:total(:percent) eta: :eta",
            total = nboot
        )
        pb_refit$tick(0)
        on.exit(pb_refit$terminate(), add = TRUE)
        bootstrap_log(verbosity, "basic", "progress bar enabled for refit loop")
    }

    reps <- if (return == "models") vector("list", nboot) else NULL
    coef_draws <- NULL

    for (i in seq_len(nboot)) {
        t_rep <- proc.time()
        bootstrap_log(verbosity, "detailed", "replicate ", i, "/", nboot, ": begin")

        intensity_sim <- intensity_sim_fixed
        outcome_sim <- outcome_sim_fixed

        if (is.null(intensity_sim) && !is.null(intensity_model)) {
            bootstrap_log(verbosity, "detailed", "replicate ", i, ": sampling intensity coefficients")
            sc <- get_bootstrap_coeffs(
                intensity_model,
                sample_coefficients = sample_coefficients,
                model_label = "intensity model"
            )
            intensity_sim <- make_parametric_intensity_simulator(intensity_model, sc)
        }

        if (is.null(outcome_sim) && !is.null(outcome_model)) {
            bootstrap_log(verbosity, "detailed", "replicate ", i, ": sampling outcome coefficients")
            sc2 <- get_bootstrap_coeffs(
                outcome_model,
                sample_coefficients = sample_coefficients,
                model_label = "outcome model"
            )
            outcome_sim <- make_parametric_single_index_simulator(outcome_model, sc2)
        }

        args <- simulate_args
        if (!is.null(intensity_sim)) {
            args$intensity_fn <- intensity_sim
            if (!is.null(intensity_bound_fixed)) {
                args$intensity_bound <- intensity_bound_fixed
            }
            if (is.null(args$intensity_bound)) {
                args$intensity_bound <- compute_safe_intensity_bound(intensity_sim, simulate_args)
            }
        }
        if (!is.null(outcome_sim)) {
            args$outcome_simulator <- outcome_sim
        }

        bootstrap_log(verbosity, "detailed", "replicate ", i, ": calling simulate_SensIAT_data with n_subjects=", args$n_subjects, ", End=", args$End)
        sim <- do.call(simulate_SensIAT_data, args)
        bootstrap_log(verbosity, "detailed", "replicate ", i, ": simulated rows=", nrow(sim))

        t_fit <- proc.time()
        bootstrap_log(verbosity, "detailed", "replicate ", i, ": fitting marginal model")
        replicated_model <-
            sim |>
            fit_SensIAT_within_group_model(
                outcome_modeler = within_group_model$outcome_modeler,
                id = !!within_group_model$variables$id,
                outcome = !!within_group_model$variables$outcome,
                time = !!within_group_model$variables$time,
                knots = unique(within_group_model$base@knots),
                alpha = within_group_model$alpha,
                End = within_group_model$End,
                intensity.args = within_group_model$args$intensity,
                outcome.args = within_group_model$args$outcome,
                influence.args = within_group_model$args$influence,
                spline.degree = within_group_model$base@order - 1L,
                add.terminal.observations = FALSE,
                link = within_group_model$link %||% "identity",
                loss = within_group_model$loss %||% "lp_mse",
                term2_method = within_group_model$term2_method %||% "fast"
            )
        bootstrap_log(
            verbosity,
            "detailed",
            "replicate ", i,
            ": fit completed in ", sprintf("%.2f", (proc.time() - t_fit)[[3]]), " s"
        )

        if (return == "models") {
            if (isTRUE(prune_models)) {
                replicated_model <- prune(replicated_model)
                bootstrap_log(verbosity, "detailed", "replicate ", i, ": model pruned")
            }
            reps[[i]] <- replicated_model
            bootstrap_log(verbosity, "detailed", "replicate ", i, ": stored fitted model")
        } else {
            if (is.null(coef_draws)) {
                coef_draws <- purrr::map(replicated_model$coefficients, function(beta) {
                    matrix(
                        NA_real_,
                        nrow = nboot,
                        ncol = length(beta),
                        dimnames = list(NULL, names(beta))
                    )
                })
                bootstrap_log(
                    verbosity,
                    "detailed",
                    "initialized coefficient storage for ", length(coef_draws), " alpha value(s)"
                )
            }
            for (j in seq_along(coef_draws)) {
                coef_draws[[j]][i, ] <- replicated_model$coefficients[[j]]
            }
            bootstrap_log(verbosity, "detailed", "replicate ", i, ": stored marginal coefficients")
        }

        if (progress && exists("pb_refit")) {
            try(pb_refit$tick(), silent = TRUE)
        }
        if (!is.null(gc_every) && gc_every > 0L && (i %% gc_every) == 0L) {
            gc(FALSE)
            bootstrap_log(verbosity, "detailed", "replicate ", i, ": explicit gc(FALSE) completed")
        }

        bootstrap_log(
            verbosity,
            "detailed",
            "replicate ", i, "/", nboot,
            ": finished (elapsed=", sprintf("%.2f", (proc.time() - t_rep)[[3]]), " s)"
        )
    }

    if (return == "models") {
        bootstrap_log(verbosity, "basic", "completed refit bootstrap; returning model list")
        return(reps)
    }

    rslt <- structure(
        list(
            bootstrap_coefficients = coef_draws,
            original_coefficients = within_group_model$coefficients,
            original_coefficient.variance = within_group_model$coefficient.variance,
            alpha = within_group_model$alpha,
            base = within_group_model$base,
            link = within_group_model$link %||% "identity",
            variables = within_group_model$variables,
            nboot = nboot,
            call = match.call(expand.dots = TRUE)
        ),
        class = "SensIAT_withingroup_bootstrap_results"
    )
    bootstrap_log(verbosity, "basic", "completed refit bootstrap; returning coefficient-only bootstrap results")
    rslt
}
