#pragma once

#include <Rcpp.h>
using namespace Rcpp;

// Forward declarations
struct IntegrationState;

/**
 * Vectorized adaptive Simpson's integration for multiple alpha values simultaneously
 * 
 * This function integrates the term2 influence function across all alpha values
 * at once, sharing the expensive computations (PMF estimation, weight functions)
 * while maintaining separate convergence criteria for each alpha.
 * 
 * @param compute_expected_values_fn R function computing E[Y*exp(alpha*Y)] and E[exp(alpha*Y)]
 * @param impute_fn R function for data imputation at time t
 * @param weight_fn R function computing weight W(t, beta)  
 * @param marginal_mean_fn R function computing marginal mean mu(t)
 * @param alpha_vec Vector of alpha sensitivity parameters
 * @param tmin Lower integration bound
 * @param tmax Upper integration bound
 * @param patient_data Patient data object passed to R functions
 * @param tol Convergence tolerance (applied per alpha)
 * @return List with integration results for each alpha
 */
List integrate_term2_vectorized_alpha(
    Function compute_expected_values_fn,
    Function impute_fn,
    Function weight_fn,
    Function marginal_mean_fn,
    NumericVector alpha_vec,
    double tmin,
    double tmax,
    SEXP patient_data,
    double tol
);

/**
 * Helper function to evaluate the term2 integrand for all alphas at time t
 * 
 * @param t Time point
 * @param compute_expected_values_fn R function for expected values
 * @param impute_fn R function for data imputation
 * @param weight_fn R function for weights
 * @param marginal_mean_fn R function for marginal means
 * @param alpha_vec Vector of alpha values
 * @param patient_data Patient data
 * @return NumericMatrix with integrand values (one column per alpha)
 */
NumericMatrix evaluate_integrand_all_alphas(
    double t,
    Function compute_expected_values_fn,
    Function impute_fn,
    Function weight_fn,
    Function marginal_mean_fn,
    NumericVector alpha_vec,
    SEXP patient_data
);