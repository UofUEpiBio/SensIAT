globalVariables(c(
    "..visit_number..", "term1", "term2", "IF", "IF_ortho",
    "estimate", "Beta_hat", "Var_beta"
))

#' Produce fitted model for group (treatment or control)
#'
#' Produces a fitted model that may be used to produce estimates of mean and
#' variance for the given group.
#'
#' This function should be agnostic to whether it is being provided a
#' treatment or control group.
#'
#' @param group.data The data for the group that is being analyzed.
#'          Preferably passed in as a single `tibble` that internally is
#'          subsetted/filtered as needed.
#' @param outcome_modeler function for fitting the outcome model.
#'          Called with a formula, data argument and `outcome.args` list.
#' @param knots knot locations for defining the spline basis.
#' @param id The variable that identifies the patient.
#' @param outcome The variable that contains the outcome.
#' @param time The variable that contains the time.
#' @param alpha The sensitivity parameter.
#' @param End The end time for this data analysis, we need to set the default value as the
#'          max value of the time.
#' @param intensity.args A list of optional arguments for intensity model.
#'          See the Intensity Arguments section.
#' @param outcome.args parameters as needed passed into the `outcome_modeler`.
#'          One special element may be `'model.modifications'` which, if present,
#'          should be a formula that will be used to modify the outcome model per, [update.formula].
#' @param influence.args A list of optional arguments used when computing the influence.
#'         See the Influence Arguments section.
#' @param spline.degree The degree of the spline basis.
#' @param add.terminal.observations Logical indicating whether to add terminal
#'          observations to the data.  If TRUE, data may not contain any `NA`s.
#'          if FALSE, data will be assumed to already include the terminal
#'          observations
#' @param loss The loss function to use. Options are `"lp_mse"` (default) and `"quasi-likelihood"`.
#'          Only used when `link` is not `"identity"`.
#' @param link The link function to use. Options are `"identity"` (default), `"log"`, and `"logit"`.
#'          When `"identity"`, the original analytic formula is used; otherwise, the generalized
#'          iterative solver is used.
#' @param term2_method Method for computing term2 influence components. Options are:
#'          `"fast"` (default), `"original"`, `"fixed_grid"`, `"seeded_adaptive"`, `"gauss_legendre"`.
#'          Only used when `link` is not `"identity"`.
#' @param impute_data A function that takes `(t, df)` and returns the imputed data at time `t`.
#'          If `NULL` (default), a standard imputation function is generated automatically.
#'
#' @return
#'  A `SensIAT_within_group_model` object containing:
#'  \itemize{
#'    \item `models`: List with `intensity` and `outcome` sub-models
#'    \item `data`: The transformed data used for fitting
#'    \item `variables`: List of variable names (id, outcome, time)
#'    \item `End`: The end time for the analysis
#'    \item `influence`: Influence function values by patient
#'    \item `alpha`: Sensitivity parameter values
#'    \item `coefficients`: Estimated spline coefficients for each alpha
#'    \item `coefficient.variance`: Variance of coefficients for each alpha
#'    \item `base`: The [`SplineBasis`][orthogonalsplinebasis::SplineBasis] object
#'    \item `V_inverse`: Inverse of the Gram matrix
#'    \item `link`: The link function used
#'    \item `loss`: The loss function used
#'    \item `term2_method`: The term2 integration method used
#'  }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Example 1: Default identity link (original analytic solution)
#' model <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         End = 830,
#'         knots = c(60, 260, 460),
#'     )
#'
#' # Example 2: Log link with generalized iterative solver
#' model_log <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
#'         alpha = 0,
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         End = 830,
#'         knots = c(60, 260, 460),
#'         link = "log",
#'         loss = "lp_mse"
#'     )
#' }
fit_SensIAT_within_group_model <- function(group.data,
                                           outcome_modeler,
                                           id,
                                           outcome,
                                           time,
                                           knots,
                                           alpha = 0,
                                           End = NULL,
                                           intensity.args = list(),
                                           outcome.args = list(),
                                           influence.args = list(),
                                           spline.degree = 3,
                                           add.terminal.observations = TRUE,
                                           loss = c("lp_mse", "quasi-likelihood"),
                                           link = c("identity", "log", "logit"),
                                           term2_method = c("fast", "original", "fixed_grid", "seeded_adaptive", "gauss_legendre"),
                                           impute_data = NULL) {
    ###### Input clean and capture -------------------------------------------
    # Variables
    id.var <- ensym(id)
    outcome.var <- ensym(outcome)
    time.var <- ensym(time)
    vars <- list(
        outcome = outcome.var,
        id = id.var,
        time = time.var
    )


    # Pass through Argument Lists
    intensity.args <- match.names(intensity.args, c("model.modifications", "bandwidth"), FALSE)
    outcome.args <- match.names(outcome.args, c("model.modifications"))
    influence.args <- match.names(influence.args, c("tolerance"))

    # Function
    outcome_modeler <- match.fun(outcome_modeler)
    if (is.null(End)) {
        End <- rlang::expr(max({{ time }}, na.rm = TRUE) + 1) %>%
            rlang::eval_tidy(data = group.data, env = parent.frame())
    }

    # Match new generalized parameters
    loss <- match.arg(loss)
    link <- match.arg(link)
    term2_method <- match.arg(term2_method)
    
    # Determine whether to use generalized solver
    use_generalized <- link != "identity"

    ###### Create Model Frame -----------------------------------------------
    outcome.extra.vars <- all.vars(outcome.args$model.modifications) |>
        setdiff(c("..id..", "..time..", "..prev_outcome..", "..delta_time..", "..visit_number..", "..outcome.."))
    intensity.extra.vars <- all.vars(intensity.args$model.modifications) |>
        setdiff(c("..id..", "..time..", "..prev_outcome..", "..delta_time..", "..visit_number..", "..outcome.."))

    data_all_with_transforms <- group.data |>
        select(
            !!id.var, !!time.var, !!outcome.var,
            any_of(outcome.extra.vars),
            any_of(intensity.extra.vars)
        ) |>
        prepare_SensIAT_data(
            id.var = !!id.var,
            time.var = !!time.var,
            outcome.var = !!outcome.var,
            End = End,
            add.terminal.observations = add.terminal.observations
        )

    ######   Intensity model  ##################################################
    #' @section Intensity Arguments:
    #' The `intensity.args` list may contain the following elements:
    #'
    #' * **`model.modifications`** A formula that will be used to modify the intensity model from it's default, per [update.formula].
    #' * **`kernel`** The kernel function for the intensity model. Default is the Epanechnikov kernel.
    #' * **`bandwidth`** The bandwidth for the intensity model kernel.
    #'
    followup_data <- data_all_with_transforms |>
        filter(.data$..time.. > 0, !is.na(.data$..prev_outcome..))

    intensity.formula <- Surv(..prev_time.., ..time.., !is.na(..outcome..)) ~ ..prev_outcome.. + strata(..visit_number..)
    if (!is.null(intensity.args$model.modifications)) {
        intensity.formula <- update.formula(intensity.formula, intensity.args$model.modifications)
    }

    intensity.model <- rlang::inject(
        coxph(intensity.formula,
            id = ..id.., data = followup_data,
            !!!purrr::discard_at(intensity.args, c("bandwidth", "kernel", "model.modifications"))
        )
    )
    # baseline_intensity_all <-
    #     estimate_baseline_intensity(
    #         intensity.model = intensity.model,
    #         data = followup_data,
    #         kernel = intensity.args$kernel %||% \(x) 0.75*(1 - (x)**2) * (abs(x) < 1),
    #         bandwidth = intensity.args$bandwidth
    #     )
    attr(intensity.model, "bandwidth") <- intensity.args$bandwidth
    attr(intensity.model, "kernel") <- intensity.args$kernel %||% \(x) 0.75 * (1 - (x)**2) * (abs(x) < 1)
    # followup_data$baseline_intensity <- baseline_intensity_all$baseline_intensity

    ######   Outcome model #####################################################
    outcome.formula <-
        ..outcome..~ -1 +
        ns(..prev_outcome.., df = 3) +
        scale(..time..) +
        scale(..delta_time..)
    if (!is.null(outcome.args$model.modifications)) {
        outcome.formula <- update.formula(outcome.formula, outcome.args$model.modifications)
    }

    outcome.model <- rlang::inject(
        outcome_modeler(outcome.formula,
            data = filter(followup_data, !is.na(..outcome..)),
            !!!purrr::discard_at(outcome.args, "model.modifications")
        )
    )


    # Compute value of the influence function: -----------------------------
    #' @section Influence Arguments:
    #'
    #' The `influence.args` list may contain the following elements:
    #'
    #' * **`method`** The method for integrating, adaptive or fixed quadrature. Default is `'adaptive'`.
    #' * **`tolerance`** The tolerance when using adaptive quadrature.
    #' * **`delta`** The bin width for fixed quadrature.
    #' * **`resolution`** alternative to `delta` by specifying the number of bins.
    #' * **`fix_discontinuity`** Whether to account for the discontinuity in the influence at observation times.
    
    # Route to either the original (identity link) or generalized marginal model
    if (use_generalized) {
        # Create default impute_data function if not provided
        if (is.null(impute_data)) {
            # Capture outcome variable name from formula
            outcome_var_name <- as.character(rlang::f_lhs(formula(outcome.model)))
            
            impute_data <- function(t, df) {
                # Set up the data with appropriate lag variables
                data_wl <- df |>
                    dplyr::mutate(
                        ..prev_time.. = .data$..time..,
                        ..prev_outcome.. = .data[[outcome_var_name]],
                        ..delta_time.. = 0
                    )
                extrapolate_from_last_observation(t, data_wl, "..time..", slopes = c("..delta_time.." = 1))
            }
        }
        
        marginal_model <- rlang::inject(fit_SensIAT_marginal_mean_model_generalized(
            data = data_all_with_transforms,
            time = data_all_with_transforms$..time..,
            id = data_all_with_transforms$..id..,
            alpha = alpha,
            knots = knots,
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            impute_data = impute_data,
            loss = loss,
            link = link,
            spline.degree = spline.degree,
            term2_method = term2_method,
            !!!influence.args
        ))
        
        # Normalize return structure (generalized uses V.inverse, original uses V_inverse)
        marginal_model$V_inverse <- marginal_model$V.inverse %||% marginal_model$V_inverse
        marginal_model$coefficient.variance <- marginal_model$coefficient.variance %||% marginal_model$coefficient_variance
    } else {
        marginal_model <- rlang::inject(fit_SensIAT_marginal_mean_model(
            data_all_with_transforms,
            alpha = alpha,
            knots = knots,
            intensity.model = intensity.model,
            outcome.model = outcome.model,
            spline.degree = spline.degree,
            !!!influence.args,
            time.vars = c("..delta_time..")
        ))
    }

    structure(
        list(
            models = list(
                intensity = intensity.model,
                outcome = outcome.model
            ),
            data = data_all_with_transforms,
            variables = vars,
            End = End,
            influence = marginal_model$influence,
            outcome_modeler = outcome_modeler,
            alpha = alpha,
            coefficients = marginal_model$coefficients,
            coefficient.variance = marginal_model$coefficient.variance,
            args = list(
                intensity = intensity.args,
                outcome = outcome.args,
                influence = influence.args
            ),
            base = marginal_model$base,
            V_inverse = marginal_model$V_inverse,
            link = link,
            loss = loss,
            term2_method = term2_method
        ),
        class = "SensIAT_within_group_model",
        call = match.call(expand.dots = TRUE)
    )
}
