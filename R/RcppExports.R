# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Runs a *basic* implementation of the "NW" function with the "K2_Biweight" kernel, just as a
#' proof-of-concept.
#'
#' @param Xb    a vector (expected to be about 500 elements)
#' @param Y     a vector (same size as `Xb`)
#' @param xb    a vector
#' @param y_seq a vector
#' @param h     a scalar, the bandwidth of kernel
#'
#' @return A matrix of the same size as `xb` by `y_seq`.
#'
#' @keywords internal
pcoriaccel_NW_basic <- function(Xb, Y, xb, y_seq, h) {
    .Call(`_SensIAT_pcoriaccel_NW_basic`, Xb, Y, xb, y_seq, h)
}

#' Runs an optimized implementation of the "NW" function.
#'
#' @param Xb    a vector (expected to be about 500 elements)
#' @param Y     a vector (same size as `Xb`)
#' @param xb    a vector
#' @param y_seq a vector
#' @param h     a scalar, the bandwidth of kernel
#' @param kernel a string, denoting the kernel function to use, either `"dnorm"`, `"K2_Biweight"`, or `"K4_Biweight"`
#'
#' @return A matrix of the same size as `xb` by `y_seq`.
#'
#' @keywords internal
pcoriaccel_NW <- function(Xb, Y, xb, y_seq, h, kernel = "K2_Biweight") {
    .Call(`_SensIAT_pcoriaccel_NW`, Xb, Y, xb, y_seq, h, kernel)
}

#' Multiplies two matrices.  If the first argument is a vector, it is interpreted as a row vector.
#' Otherwise, if the second argument is a vector, it is interpreted as a column vector.
#' @param matrA first matrix
#' @param matrB second matrix
#' @return The product of `matrA` and `matrB`.
#' @keywords internal
#' @noRd
pcoriaccel_mmul <- function(matrA, matrB) {
    .Call(`_SensIAT_pcoriaccel_mmul`, matrA, matrB)
}

#' Inner product (dot product) of two vectors.
#' @param vecA first vector
#' @param vecB second vector
#' @return scalar product of `vecA` and `vecB`
#' @keywords internal
#' @noRd
pcoriaccel_inner_prod <- function(vecA, vecB) {
    .Call(`_SensIAT_pcoriaccel_inner_prod`, vecA, vecB)
}

#' Outer sum of two vectors.
#' @param vecA first vector
#' @param vecB second vector
#' @return Matrix where each element of `vecA`(row) is added to the element of `vecB`(column).
#' @keywords internal
#' @noRd
pcoriaccel_outer_sum <- function(vecA, vecB) {
    .Call(`_SensIAT_pcoriaccel_outer_sum`, vecA, vecB)
}

#' Outer product of two vectors.
#' @param vecA first vector
#' @param vecB second vector
#' @return matrix of the outer product of vectors `vecA` and `vecB`.
#' @keywords internal
#' @noRd
pcoriaccel_outer_prod <- function(vecA, vecB) {
    .Call(`_SensIAT_pcoriaccel_outer_prod`, vecA, vecB)
}

#' Returns the unique elements of a vector, sorted in ascending order.
#' @param vec the vector
#' @return `sort(unique(vec))`
#' @keywords internal
#' @noRd
pcoriaccel_sorted_unique <- function(vec) {
    .Call(`_SensIAT_pcoriaccel_sorted_unique`, vec)
}

#' Runs an optimized implementation of the `compute_influence_term_2_quadv_sim_via_matrix`
#' function.
#'
#' @param X              Matrix of all covariates, transformed as necessary by model
#' @param Y              Vector of all outcomes (same length as a column of `X`)
#' @param times          Vector of observation times for individual
#' @param individual_X   Matrix of covariates for individual rows correspond to times prepared for
#'                       inferences for integration.
#' @param x_slope        Vector of numeric(length(beta)) indicating how
#' @param alpha          Vector of sensitivity parameters
#' @param beta           Vector of coefficients of the outcome model
#' @param spline_basis   Spline basis object (`orthogonalsplinebasis::SplineBasis`)
#' @param bandwidth      Bandwidth for the kernel density estimate of the outcome model.
#' @param tol            Tolerance for integration
#' @param kernel         Kernel function to use for the kernel density estimate
#'
#' @return integration result
#'
#' @keywords internal
pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix <- function(X, Y, times, individual_X, x_slope, alpha, beta, spline_basis, bandwidth, tol = 0.0001220703, kernel = "K2_Biweight") {
    .Call(`_SensIAT_pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix`, X, Y, times, individual_X, x_slope, alpha, beta, spline_basis, bandwidth, tol, kernel)
}

#' Directly estimate the probability mass function of Y.
#'
#' @param Xb Numeric vector of individual linear predictors from the data
#' @param Y Numeric vector of individual responses from the data
#' @param xi value of the individuals linear predictor at the point of estimation
#' @param y_seq Numeric vector of unique values of `Y`.
#' @param h bandwidth of the kernel
#' @param kernel character string specifying the kernel to use, either `"dnorm"`, `"K2_Biweight"`, or `"K4_Biweight"`
#'
pcoriaccel_estimate_pmf <- function(Xb, Y, xi, y_seq, h, kernel = "K2_Biweight") {
    .Call(`_SensIAT_pcoriaccel_estimate_pmf`, Xb, Y, xi, y_seq, h, kernel)
}

#' Integrate function using adaptive Simpson quadrature.
#'
#' @param integrand   The integrand, must take scalar argument, may return scalar, vector, or matrix.
#' @param lo          Lower integration bound
#' @param hi          Upper integration bound
#' @param tol         Tolerance for integration, default `.Machine$double.eps^(1/2)`
#'
#' @return integration result, list with elements `$Q` (the integral estimate), `$fcnt` (the number
#' of function evaluations), and `$estim.prec` (a (pessimistic) estimate of the precision).
#'
#' @keywords internal
pcoriaccel_integrate_simp <- function(integrand, lo, hi, tol = 1.490116e-08) {
    .Call(`_SensIAT_pcoriaccel_integrate_simp`, integrand, lo, hi, tol)
}

#' Compiled version of `evaluate_basis()` function
#'
#' @param spline_basis   The spline basis, S4 class `orthogonalsplinebasis::SplineBasis`
#' @param x              The point to evaluate
#'
#' @return Vector of the basis functions evaluated at x.
#'
pcoriaccel_evaluate_basis <- function(spline_basis, x) {
    .Call(`_SensIAT_pcoriaccel_evaluate_basis`, spline_basis, x)
}

