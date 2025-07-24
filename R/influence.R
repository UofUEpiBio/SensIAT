
#' Compute Influence Terms
#'
#' This function computes the influence terms for the marginal outcome model sensitivity analysis.
#' It is a generic function that can handle different types of outcome models.
#'
#' @param outcome.model The outcome model fitted to the data.
#' @param intensity.model The intensity model fitted to the data.
#' @param alpha A numeric vector representing the sensitivity parameter.
#' @param data A data frame containing the observations.
#' @param id A variable representing the patient identifier.
#' @param base A spline basis object.
#' @param ... Additional arguments passed to the method.
#'
#' @export
compute_influence_terms <-
function(
    outcome.model, # outcome model
    intensity.model, # The intensity model
    alpha, # Sensitivity, singular alpha value
    data, #< vector of times for all observations
    ...
){
    UseMethod("compute_influence_terms", outcome.model)
}

#' @describeIn compute_influence_terms Generic method which covers lm and glm outcome models.
#' @export
compute_influence_terms.default <-
function(outcome.model, intensity.model, alpha, data, id, base, ...){
    rlang::abort("compute_influence_terms is not implemented for this outcome model type. Please use a method that supports your model type.")
}
