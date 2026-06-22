## Parametric bootstrap orchestration and helpers

#' Sample parametric coefficients from a fitted model using asymptotic MVN
#'
#' @param model A fitted model with `coef()` and `vcov()` methods
#' @return Named numeric vector of sampled coefficients
sample_parametric_coeffs <- function(model) {
    coefs <- tryCatch(coef(model), error = function(e) NULL)
    if (is.null(coefs)) stop("Model does not provide coef()")
    vc <- tryCatch(vcov(model), error = function(e) NULL)
    if (is.null(vc)) stop("Model does not provide vcov(); cannot sample coefficients")
    # Use MASS::mvrnorm
    samp <- as.numeric(MASS::mvrnorm(1, mu = as.numeric(coefs), Sigma = vc))
    names(samp) <- names(coefs)
    return(samp)
}


#' Build a parametric intensity simulator using sampled coefficients
#'
#' Supports `coxph` objects: sampled coefficients are used to compute linear predictor
#' while baseline hazard is taken from `survival::basehaz()`.
#'
make_parametric_intensity_simulator <- function(intensity_model, sampled_coef, covariate_mapping = NULL) {
    if (!inherits(intensity_model, "coxph")) stop("Only coxph intensity_model supported for parametric sampling currently")

    model_terms <- delete.response(terms(intensity_model))
    bh <- tryCatch(survival::basehaz(intensity_model, centered = FALSE), error = function(e) NULL)
    if (is.null(bh) || nrow(bh) == 0) {
        warning("Could not extract baseline hazard from coxph; using constant hazard 0.01")
        return(function(t, prev_outcome, visit_num) 0.01)
    }
    times <- c(0, bh$time)
    cumhaz <- c(0, if (!is.null(bh$hazard)) bh$hazard else bh$haz)
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
            coef_names <- names(sampled_coef)
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
make_parametric_single_index_simulator <- function(outcome_model, sampled_coef, covariate_mapping = NULL) {
    if (!inherits(outcome_model, "SensIAT::Single-index-outcome-model")) stop("outcome_model must be a fitted single-index outcome model")
    bandwidth <- if (!is.null(outcome_model$bandwidth)) outcome_model$bandwidth else stop("bandwidth missing in outcome_model")
    data_orig <- outcome_model$data
    default_names <- c(prev_outcome = "prev_outcome", time = "time", delta_time = "delta_time")
    if (!is.null(covariate_mapping)) {
        if (!is.character(covariate_mapping) || is.null(names(covariate_mapping))) stop("covariate_mapping must be a named character vector")
        default_names[names(covariate_mapping)] <- covariate_mapping
    }

    model_terms <- delete.response(terms(outcome_model))
    y_resp <- model.response(model.frame(outcome_model))
    y_seq <- sort(unique(y_resp))

    function(prev_outcome, time, delta_time, newdata = NULL) {
        if (!is.null(newdata)) {
            if (!is.data.frame(newdata) || nrow(newdata) < 1) stop("newdata must be a data.frame with at least one row")
            nd <- newdata[1, , drop = FALSE]
        } else {
            nd <- data.frame(prev_outcome = prev_outcome, time = time, delta_time = delta_time, stringsAsFactors = FALSE)
            names(nd) <- unname(default_names)
        }

        X_new <- tryCatch(model.matrix(model_terms, data = nd), error = function(e) stop("Failed to build model matrix for outcome newdata: ", conditionMessage(e)))
        X_orig <- model.matrix(model_terms, data = data_orig)
        lp0 <- as.vector(X_orig %*% as.numeric(sampled_coef[colnames(X_orig)]))
        xb <- as.vector(X_new %*% as.numeric(sampled_coef[colnames(X_new)]))

        Fhat <- pcoriaccel_NW(Xb = lp0, Y = y_resp, xb = xb[1], y_seq = y_seq, h = bandwidth, kernel = attr(outcome_model, "kernel"))
        Fhat_vec <- as.vector(Fhat)
        pmf <- pmax(0, diff(c(0, Fhat_vec)))
        if (sum(pmf) <= 0) return(mean(y_resp, na.rm = TRUE))
        pmf <- pmf / sum(pmf)
        sampled <- sample(y_seq, size = 1, prob = pmf)
        return(sampled)
    }
}


#' Parametric bootstrap orchestration
#'
#' @param nboot Number of bootstrap replicates
#' @param intensity_model Fitted intensity model (coxph) or function
#' @param outcome_model Fitted outcome model (single-index) or NULL
#' @param simulate_args List of arguments to pass to `simulate_SensIAT_data()` (e.g., n_subjects, End, intensity_bound)
#' @param seed Optional seed for reproducibility
#' @return A list of simulated datasets (length nboot)
parametric_bootstrap <- function(nboot = 100,
                                 intensity_model = NULL,
                                 outcome_model = NULL,
                                 simulate_args = list(),
                                 seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    results <- vector("list", nboot)

    for (b in seq_len(nboot)) {
        # Sample coefficients and build simulators as available
        intensity_sim <- NULL
        outcome_sim <- NULL

        if (!is.null(intensity_model)) {
            sc <- sample_parametric_coeffs(intensity_model)
            intensity_sim <- make_parametric_intensity_simulator(intensity_model, sc)
        }

        if (!is.null(outcome_model)) {
            sc2 <- sample_parametric_coeffs(outcome_model)
            outcome_sim <- make_parametric_single_index_simulator(outcome_model, sc2)
        }

        # Merge supplied simulate_args and inject simulators
        args <- simulate_args
        if (!is.null(intensity_sim)) {
            args$intensity_fn <- intensity_sim
            # require an intensity_bound; try to set from simulate_args or compute a safe bound
            if (is.null(args$intensity_bound)) {
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
            }
        }
        if (!is.null(outcome_sim)) args$outcome_simulator <- outcome_sim

        # Call simulation
        sim <- do.call(simulate_SensIAT_data, args)
        results[[b]] <- sim
    }
    return(results)
}

#' Parametric bootstrap for a within-group SensIAT model
#'
#' @param nboot Number of bootstrap replicates.
#' @param within_group_model A fitted `SensIAT_within_group_model` object.
#' @param simulate_args List of arguments to pass to `simulate_SensIAT_data()`.
#'        If not specified, `End`, `n_subjects`, `initial_outcome_mean`, and
#'        `initial_outcome_sd` are inferred from the fitted model.
#' @param seed Optional seed for reproducibility.
#' @return A list of simulated datasets, one per bootstrap replicate.
#' @export
parametric_bootstrap_within_group <- function(nboot = 100,
                                              within_group_model,
                                              simulate_args = list(),
                                              seed = NULL) {
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

    parametric_bootstrap(
        nboot = nboot,
        intensity_model = intensity_model,
        outcome_model = outcome_model,
        simulate_args = simulate_args,
        seed = seed
    )
}
