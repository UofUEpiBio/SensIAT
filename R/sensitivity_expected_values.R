

#' Compute Conditional Expected Values based on Outcome Model
#'
#' @param model An object representing the output of the outcome model.
#' @param alpha The sensitivity parameter
#' @param new.data Data to compute conditional means for, defaults to the model frame for the fitted model.
#' @param ... passed onto methods.
#'
#' @details
#' Compute the conditional expectations needed for predictions in the models.
#' Two additional values/expectations are computed:
#'
#' * `$E \big[ Y(t) \exp \{  \alpha Y(t) \}   | A(t)=1, \bar{O}(t) \big]$`, returned as `E_Yexp_alphaY`, and
#' * `$E \big[ \exp \{ \alpha Y(t) \} \  | A(t)=1, \bar{O}(t) \big]$`, returned as `E_exp_alphaY`.
#'
#'
#' @return The `new.data` frame with additional columns `E_Yexp_alphaY`, and `E_exp_alphaY` appended.
#' @export
sensitivity_expected_values <-
function(
    model,
    alpha = 0, #< sensitivity parameter
    new.data = model.frame(model),
    ...){
    UseMethod('pcori_conditional_means')
}
