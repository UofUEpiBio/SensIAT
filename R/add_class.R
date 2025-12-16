#' Adds an S3 class to an object
#'
#' This is a utility function for working with the S3 class system only.
#' It does not work with S4 or R6 class systems.
#'
#' @param x An object to which the class should be added.
#' @param class A character vector of class names to be added.
#'
#' @return The object with the new class(es) prepended to existing classes.
#'
#' @keywords internal
#' @examples
#' # Internal use only
#' obj <- add_class(list(a = 1), "my_class")
#' class(obj)  # c("my_class", "list")
add_class <- function(x, class) {
    if (!is.character(class)) {
        stop("class must be a character vector")
    }

    # Check if the object already has the class
    existing_classes <- class(x)

    # If the class is already present, return the object as is
    if (class %in% existing_classes) {
        return(x)
    }

    # Add the new class to the object
    class(x) <- c(class, existing_classes)

    return(x)
}
