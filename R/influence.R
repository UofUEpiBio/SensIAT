
#' Compute Influence Terms
#'
#' @export
compute_influence_terms <-
function(
    data_all_with_transforms, #< vector of times for all observations
    base, # Spline basis
    alpha, # Sensitivity, singular alpha value
    outcome.model, # outcome model
    intensity.model, # The intensity model
    ...
){
    UseMethod("compute_influence_terms", outcome.model)
}

#' @describeIn compute_influence_terms Generic method which covers lm and glm outcome models.
#' @export
compute_influence_terms.default <-
function(
    time, #< vector of times for all observations
    base, # Spline basis
    alpha, # Sensitivity, recurses if length > 1
    outcome.model, # outcome model
    intensity.model, # The intensity model
    ...
){

    term_1 <- compute_influence_term_1(
        time, base,
        outcome.model = outcome.model,
        intensity.model = intensity.model,
        alpha = alpha
    )

    term_2 <- compute_influence_term_2
}
compute_influence_term_1 <-
    function(
        time, #< vector of times for all observations
        base, # Spline basis
        outcome.model, # outcome model
        intensity.model, # The intensity model
        alpha, # Sensitivity, singular alpha value
        ...
    ){
        term_1_scalar <- compute_influence_term_1_scalar(
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            alpha = alpha
        )
        pcoriaccel_evaluate_basis_mat(base, time) * c(term_1_scalar)
    }


compute_influence_term_1_scalar <-
    function(
        outcome.model, # outcome model
        intensity.model, # The intensity model
        alpha, # Sensitivity, singular alpha value
        ...
    ){
        expected.values <- sensitivity_expected_values(outcome.model, alpha)

        outcome = model.response(model.frame(outcome.model))

        rho <- compute_influence_term_1_rho(
            outcome,
            baseline_intensity =
                estimate_baseline_intensity(intensity.model)$baseline_intensity,
            predict(intensity.model, type='lp', reference='zero'),
            alpha,
            expected.values$E_exp_alphaY
        )

        (outcome - expected.values$E_Yexp_alphaY/expected.values$E_exp_alphaY) /
            rho
    }

compute_influence_term_1_rho <-
    function(
        outcome,
        baseline_intensity,
        exp_intensity_linear_predictor,
        alpha,
        E_exp_alphaY
    ){
        baseline_intensity *
        exp_intensity_linear_predictor *
        exp(-alpha*outcome) *
        E_exp_alphaY
    }

compute_influence_term_2 <-
    function(
        base, # Spline basis
        outcome.model,
        alpha
    ){
        lower <- base@knots[base@order]
        upper <- base@knots[length(base@knots) - base@order+1]





        pracma::quadv(\(t){
            evaluate(base, t) *
            EV <- sensitivity_expected_values(outcome.model, alpha, t)

        }, lower, upper)



    }

compute_influence_term_2_for_one_patient <-
    function(
        base, # Spline basis
        outcome.model,
        alpha
    ){

        expected_value <- function(time, lhs=FALSE){
            if(lhs){
                period <- max(which(times - time <= 0))
            } else {
                period <- max(which(times - time <  0))
            }
            xi <- individual_X[period, ,drop=FALSE] + (time - times[period]) * x_slope


            pmf <- pcoriaccel_estimate_pmf(
                Xb=Xb, Y = Y,
                xi = xi %*% beta,
                y_seq = uY,
                h= bandwidth, kernel = kernel)
            if(all(pmf ==0)) return(0)
            E_exp_alphaY <- sum( exp(alpha*uY)*pmf )

            E_Yexp_alphaY <- sum( uY*exp(alpha*uY)*pmf )

            return(E_Yexp_alphaY/E_exp_alphaY)
        }


    }
