

#' Compute Conditional Means
#'
#' @param model An object of class `SensIAT::outcome-model`
#' @param alpha Sensitivity parameter
#' @param new.data Data to compute conditional means for, defaults to the model frame for the fitted model.
#' @param ... passed onto methods.
#'
#' @details
#' Compute the conditional expectations needed for predictions in the models.
#' Three additional values/expectations are computed:
#'
#' * `$E \big[ Y(t) \exp \{  \alpha Y(t) \}   | A(t)=1, \bar{O}(t) \big]$`, returned as `E_y_past`, and
#' * `$E \big[ \exp \{ \alpha Y(t) \} \  | A(t)=1, \bar{O}(t) \big]$`, returned as `E_exp_alphaY`.
#'
#'
#' @return The `new.data` frame with additional columns `E_Y_past`, and `E_exp_alphaY` appended.
pcori_conditional_means <-
function(
    model,
    alpha = 0, #< sensitivity parameter
    new.data = model.frame(model),
    ...){
    UseMethod('pcori_conditional_means')
}

