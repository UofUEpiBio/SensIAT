#' Adds an S3 class to an object
#'
#' @param x An object to which the class should be added.
#' @param class A character vector of class names to be added.
#'
#' @export
#' @examples
#' add_class(TRUE, 'flag')
#'
#' @keywords internal
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
