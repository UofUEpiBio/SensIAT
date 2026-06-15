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
    times <- bh$time
    cumhaz <- if (!is.null(bh$hazard)) bh$hazard else bh$haz
    H0 <- stats::stepfun(times, c(0, cumhaz), right = TRUE)

    function(t, prev_outcome, visit_num) {
        nd <- data.frame(prev_outcome = prev_outcome, visit_num = visit_num)
        # Build model matrix and compute lp using sampled_coef
        X_new <- tryCatch(model.matrix(model_terms, data = nd), error = function(e) NULL)
        if (is.null(X_new)) {
            lp <- 0
        } else {
            # Align coefficients
            coef_names <- names(sampled_coef)
            cols <- colnames(X_new)
            # Some models include an intercept column '(Intercept)'
            coef_vec <- sampled_coef
            missing <- setdiff(cols, names(coef_vec))
            if (length(missing) > 0) {
                # set missing coef to 0
                coef_vec[missing] <- 0
            }
            lp <- as.numeric(X_new %*% coef_vec[cols])
        }
        dt <- max(1e-5, 1e-4 * max(1, abs(t)))
        h <- (H0(t + dt) - H0(t)) / dt
        if (!is.finite(h) || h < 0) h <- 0
        val <- h * exp(lp)
        if (!is.finite(val)) val <- 0
        return(as.numeric(val))
    }
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
            # require an intensity_bound; try to set from simulate_args or default
            if (is.null(args$intensity_bound)) args$intensity_bound <- 0.05
        }
        if (!is.null(outcome_sim)) args$outcome_simulator <- outcome_sim

        # Call simulation
        sim <- do.call(simulate_SensIAT_data, args)
        results[[b]] <- sim
    }
    return(results)
}
