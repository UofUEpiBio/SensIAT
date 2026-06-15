## Factory to create simulators from fitted Single-index outcome models
##
#' Create a simulator function from a fitted Single-index outcome model
#'
#' This factory builds a closure that can be used to simulate outcome values
#' from a fitted `SensIAT::Single-index-outcome-model`. The returned function
#' samples from the estimated conditional distribution of the outcome given
#' covariates using the Nadaraya-Watson estimator implemented in
#' `pcoriaccel_NW`.
#'
#' @param outcome_model A fitted object of class `SensIAT::Single-index-outcome-model`.
#' @param covariate_mapping Optional named character vector mapping expected
#'   covariate names (e.g. `prev_outcome`, `time`, `delta_time`) to the names
#'   used in the original model formula. If `NULL`, the factory will attempt
#'   to use `prev_outcome`, `time`, and `delta_time` directly.
#' @return A function with signature `function(prev_outcome, time, delta_time, newdata = NULL)`
#'   which returns a sampled outcome value consistent with the fitted model.
#' @export
make_single_index_simulator <- function(outcome_model, covariate_mapping = NULL) {
    if (is.null(outcome_model) || !inherits(outcome_model, "SensIAT::Single-index-outcome-model")) {
        stop("outcome_model must be a fitted Single-index outcome model (class 'SensIAT::Single-index-outcome-model').")
    }

    # Extract stored elements needed for NW sampling
    coef_vec <- coef(outcome_model)
    bandwidth <- if (!is.null(outcome_model$bandwidth)) outcome_model$bandwidth else stop("bandwidth missing in outcome_model")
    data_orig <- outcome_model$data

    # Identify default covariate names and apply mapping if provided
    default_names <- c(prev_outcome = "prev_outcome", time = "time", delta_time = "delta_time")
    if (!is.null(covariate_mapping)) {
        if (!is.character(covariate_mapping) || is.null(names(covariate_mapping))) {
            stop("covariate_mapping must be a named character vector")
        }
        default_names[names(covariate_mapping)] <- covariate_mapping
    }

    simulator <- function(prev_outcome, time, delta_time, newdata = NULL) {
        # Build newdata row either from provided newdata or from canonical inputs
        if (!is.null(newdata)) {
            if (!is.data.frame(newdata) || nrow(newdata) < 1) stop("newdata must be a data.frame with at least one row")
            nd <- newdata[1, , drop = FALSE]
        } else {
            nd <- data.frame(
                prev_outcome = prev_outcome,
                time = time,
                delta_time = delta_time,
                stringsAsFactors = FALSE
            )
            names(nd) <- unname(default_names)
        }

        # Try to build model matrix; fall back with a helpful error if variables are missing
        model_terms <- delete.response(terms(outcome_model))
        X_new <- tryCatch(
            model.matrix(model_terms, data = nd),
            error = function(e) {
                stop("Failed to construct model matrix for newdata. Check covariate names or provide 'newdata' with the same variables used to fit the outcome_model. Error: ", conditionMessage(e))
            }
        )

        # Extract observed lp and outcomes from original data
        X_orig <- model.matrix(model_terms, data = data_orig)
        lp0 <- as.vector(X_orig %*% coef_vec)
        Y <- model.response(model.frame(outcome_model))
        y_seq <- sort(unique(Y))

        # Compute fitted conditional CDF at xb using pcoriaccel_NW
        xb <- as.vector(X_new %*% coef_vec)
        # pcoriaccel_NW expects xb length possibly >1; call for first row only
        Fhat <- pcoriaccel_NW(Xb = lp0, Y = Y, xb = xb[1], y_seq = y_seq, h = bandwidth, kernel = attr(outcome_model, "kernel"))

        # Ensure numeric vector and convert to pmf
        Fhat_vec <- as.vector(Fhat)
        pmf <- pmax(0, diff(c(0, Fhat_vec)))
        if (sum(pmf) <= 0) {
            # fallback: return the unconditional mean of the observed outcomes
            return(mean(Y, na.rm = TRUE))
        }
        pmf <- pmf / sum(pmf)

        # Sample outcome from support y_seq using pmf
        sampled <- sample(y_seq, size = 1, prob = pmf)
        return(sampled)
    }

    # Attach metadata and return
    attr(simulator, "outcome_model") <- outcome_model
    class(simulator) <- c("SensIAT_single_index_simulator", class(simulator))
    return(simulator)
}


#' Create an intensity simulator from a fitted intensity model
#'
#' This factory builds a closure that computes an observation intensity
#' (hazard) at time t given covariates like `prev_outcome` and `visit_num`.
#' Supports user-supplied functions and `coxph` fitted models (approximate
#' baseline hazard via `survival::basehaz`).
#'
#' @param intensity_model A fitted model (e.g., `coxph`) or a function(t, prev_outcome, visit_num).
#' @param covariate_mapping Optional named character vector mapping expected covariate names to model variable names.
#' @return A function of signature `function(t, prev_outcome, visit_num)` returning non-negative numeric intensity.
#' @keywords internal
make_intensity_simulator <- function(intensity_model, covariate_mapping = NULL) {
    if (is.null(intensity_model)) stop("intensity_model must be provided")

    # If user provided a function already, validate and return
    if (is.function(intensity_model)) {
        f <- intensity_model
        return(function(t, prev_outcome, visit_num) {
            val <- f(t, prev_outcome, visit_num)
            if (!is.numeric(val) || length(val) != 1) stop("intensity function must return a single numeric value")
            if (val < 0) stop("intensity function returned negative value")
            return(val)
        })
    }

    # Support coxph-like objects by extracting baseline cumulative hazard
    if (inherits(intensity_model, "coxph")) {
        bh <- tryCatch(survival::basehaz(intensity_model, centered = FALSE), error = function(e) NULL)
        if (is.null(bh) || nrow(bh) == 0) {
            warning("Could not extract baseline hazard from coxph; falling back to constant hazard of 0.01")
            return(function(t, prev_outcome, visit_num) 0.01)
        }

        # bh should have columns 'time' and 'hazard' (cumulative)
        times <- bh$time
        cumhaz <- if (!is.null(bh$hazard)) bh$hazard else bh$haz
        if (is.null(cumhaz)) cumhaz <- bh$hazard

        H0 <- stats::stepfun(times, c(0, cumhaz), right = TRUE)

        intensity_fun <- function(t, prev_outcome, visit_num) {
            # Build a minimal newdata for linear predictor if model has terms
            nd <- data.frame(prev_outcome = prev_outcome, visit_num = visit_num)
            # Predict linear predictor; fall back to 0 on error
            lp <- tryCatch(as.numeric(stats::predict(intensity_model, newdata = nd, type = "lp")), error = function(e) 0)

            # Approximate hazard via small finite difference of cumulative hazard
            dt <- max(1e-5, 1e-4 * max(1, abs(t)))
            h <- (H0(t + dt) - H0(t)) / dt
            if (!is.finite(h) || h < 0) h <- 0
            val <- h * exp(lp)
            if (!is.finite(val)) val <- 0
            return(as.numeric(val))
        }

        return(intensity_fun)
    }

    stop("Unsupported intensity_model type. Provide a function or a coxph object.")
}
