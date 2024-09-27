prune_outcome_model <- function(outcome.model){
    list(
        coefficients = outcome.model$coefficients,
        bandwidth = outcome.model$bandwidth
    )
}
prune_intensity_model <- function(intensity.model){
    list(
        coefficients = intensity.model$coefficients,
        bandwidth = attr(intensity.model, 'bandwidth')
    )
}
prune_bootstrap_replication <- function(object){
    structure(
        list(
            coefficients = object$coefficients,
            coefficient.variance = object$coefficient.variance,
            intensity.bandwidth = object$intensity.bandwidth,
            outcome = prune_outcome_model(object$outcome.model),
            intensity = prune_intensity_model(object$intensity.model),
            alpha = object$alpha),
        class = c(
            'pcori_bootstrap_replication',
            paste('pruned-', oldClass(object)),
            oldClass(object)
        )
    )
}
