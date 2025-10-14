
#' @export
`prune.SensIAT::Single-index-outcome-model` <- function(tree, ...){
    tree$frame <- tree$data <- NULL
    attr(tree, 'terms') <- NULL
    add_class(tree, 'pruned-SensIAT::Single-index-outcome-model')
}

#' @export
prune.default <- function(tree, ...){
    tree
}

#' @export
prune.coxph <- function(tree, ...){
    tree
}


#' @export
prune.SensIAT_within_group_model <- function(tree, ...){
    tree$base <- NULL
    tree$data <- NULL
    tree$V_inverse <- NULL
    tree$influence <- NULL
    tree$args <- NULL

    tree$models <- lapply(X=tree$models, FUN=prune, ...)

    add_class(tree, 'pruned-SensIAT_within_group_model')
}
