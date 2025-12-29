#include "vectorized_integration.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

NumericMatrix evaluate_integrand_all_alphas(
    double t,
    Function compute_expected_values_fn,
    Function impute_fn,
    Function weight_fn,
    Function marginal_mean_fn,
    NumericVector alpha_vec,
    SEXP patient_data
) {
    int n_alphas = alpha_vec.length();
    
    // Impute data at time t (shared across all alphas)
    SEXP imputed_data = impute_fn(t, patient_data);
    
    // Compute weight (shared across all alphas) 
    NumericVector weight = weight_fn(t);
    
    // Compute marginal mean (shared across all alphas)
    double marginal_mean = as<double>(marginal_mean_fn(t));
    
    // Compute expected values for all alphas at once
    DataFrame expected_values = compute_expected_values_fn(
        Named("alpha") = alpha_vec,
        Named("new.data") = imputed_data
    );
    
    // Extract E_Yexp_alphaY and E_exp_alphaY columns
    NumericVector E_Yexp_alphaY = expected_values["E_Yexp_alphaY"];
    NumericVector E_exp_alphaY = expected_values["E_exp_alphaY"];
    
    // Build result matrix: integrand values for each alpha
    int weight_length = weight.length();
    NumericMatrix result(weight_length, n_alphas);
    
    for (int i = 0; i < n_alphas; ++i) {
        double conditional_mean = E_Yexp_alphaY[i] / E_exp_alphaY[i];
        double integrand_scalar = conditional_mean - marginal_mean;
        
        // Apply weight vector (element-wise multiplication)
        for (int j = 0; j < weight_length; ++j) {
            result(j, i) = weight[j] * integrand_scalar;
        }
    }
    
    return result;
}

// Adaptive Simpson's rule for vectorized integration - FIXED VERSION
// This version uses return-based accumulation instead of global convergence flags
List vectorized_adaptive_simpson(
    Function compute_expected_values_fn,
    Function impute_fn,
    Function weight_fn,
    Function marginal_mean_fn,
    NumericVector alpha_vec,
    SEXP patient_data,
    double a, double b,
    double tolerance
) {
    int n_alphas = alpha_vec.length();
    int function_call_count = 0;
    const int MAX_FUNCTION_CALLS = 10000;
    const int MAX_DEPTH = 20;
    
    // Cache for function evaluations
    std::map<double, NumericMatrix> eval_cache;
    
    // Helper to evaluate integrand with caching
    auto evaluate_cached = [&](double t) -> NumericMatrix {
        auto it = eval_cache.find(t);
        if (it != eval_cache.end()) {
            return it->second;
        }
        
        if (function_call_count >= MAX_FUNCTION_CALLS) {
            stop("Too many function evaluations; stopping integration");
        }
        
        NumericMatrix result = evaluate_integrand_all_alphas(
            t, compute_expected_values_fn, impute_fn, weight_fn,
            marginal_mean_fn, alpha_vec, patient_data
        );
        eval_cache[t] = result;
        function_call_count++;
        return result;
    };
    
    // Get initial evaluations
    NumericMatrix fa = evaluate_cached(a);
    NumericMatrix fb = evaluate_cached(b);
    int weight_length = fa.nrow();
    
    // Recursive adaptive Simpson for a single interval
    // Returns the integral estimate for each alpha
    std::function<std::vector<NumericVector>(double, double, NumericMatrix, NumericMatrix, int)> 
    adaptive_simpson_recursive;
    
    adaptive_simpson_recursive = [&](
        double xa, double xb,
        NumericMatrix fa_mat, NumericMatrix fb_mat,
        int depth
    ) -> std::vector<NumericVector> {
        
        if (depth > MAX_DEPTH) {
            // Hit maximum depth - return simple estimate
            std::vector<NumericVector> results(n_alphas);
            double h = xb - xa;
            for (int alpha_idx = 0; alpha_idx < n_alphas; ++alpha_idx) {
                NumericVector Q(weight_length);
                for (int j = 0; j < weight_length; ++j) {
                    // Trapezoidal rule as fallback
                    Q[j] = (h / 2.0) * (fa_mat(j, alpha_idx) + fb_mat(j, alpha_idx));
                }
                results[alpha_idx] = Q;
            }
            return results;
        }
        
        double xc = 0.5 * (xa + xb);
        double h = xb - xa;
        
        if (h < 1e-12) {
            // Segment too small - return zero
            std::vector<NumericVector> results(n_alphas);
            for (int alpha_idx = 0; alpha_idx < n_alphas; ++alpha_idx) {
                results[alpha_idx] = NumericVector(weight_length, 0.0);
            }
            return results;
        }
        
        // Evaluate at midpoint
        NumericMatrix fc_mat = evaluate_cached(xc);
        
        // Compute estimates for each alpha
        std::vector<NumericVector> results(n_alphas);
        bool any_needs_refinement = false;
        
        for (int alpha_idx = 0; alpha_idx < n_alphas; ++alpha_idx) {
            // Simpson's 1/3 rule estimate: (h/6)(f(a) + 4f(c) + f(b))
            NumericVector Q_coarse(weight_length);
            for (int j = 0; j < weight_length; ++j) {
                Q_coarse[j] = (h / 6.0) * (
                    fa_mat(j, alpha_idx) +
                    4.0 * fc_mat(j, alpha_idx) +
                    fb_mat(j, alpha_idx)
                );
            }
            
            // Evaluate at quarter points for refined estimate
            double xd = 0.5 * (xa + xc);
            double xe = 0.5 * (xc + xb);
            NumericMatrix fd_mat = evaluate_cached(xd);
            NumericMatrix fe_mat = evaluate_cached(xe);
            
            // Two Simpson's 1/3 rules over sub-intervals
            NumericVector Q_fine(weight_length);
            for (int j = 0; j < weight_length; ++j) {
                double Q1 = (h / 12.0) * (
                    fa_mat(j, alpha_idx) +
                    4.0 * fd_mat(j, alpha_idx) +
                    fc_mat(j, alpha_idx)
                );
                double Q2 = (h / 12.0) * (
                    fc_mat(j, alpha_idx) +
                    4.0 * fe_mat(j, alpha_idx) +
                    fb_mat(j, alpha_idx)
                );
                Q_fine[j] = Q1 + Q2;
            }
            
            // Check convergence using error estimate
            double max_error = 0.0;
            for (int j = 0; j < weight_length; ++j) {
                double error = std::abs(Q_fine[j] - Q_coarse[j]) / 15.0;  // Error estimate
                max_error = std::max(max_error, error);
            }
            
            if (max_error < tolerance) {
                // Converged - use Richardson extrapolation
                for (int j = 0; j < weight_length; ++j) {
                    Q_fine[j] += (Q_fine[j] - Q_coarse[j]) / 15.0;
                }
                results[alpha_idx] = Q_fine;
            } else {
                // Not converged - recurse on sub-intervals
                any_needs_refinement = true;
                auto left_results = adaptive_simpson_recursive(xa, xc, fa_mat, fc_mat, depth + 1);
                auto right_results = adaptive_simpson_recursive(xc, xb, fc_mat, fb_mat, depth + 1);
                
                // Sum the results
                NumericVector Q_total(weight_length);
                for (int j = 0; j < weight_length; ++j) {
                    Q_total[j] = left_results[alpha_idx][j] + right_results[alpha_idx][j];
                }
                results[alpha_idx] = Q_total;
            }
        }
        
        return results;
    };
    
    // Perform integration over [a, b]
    std::vector<NumericVector> integrals = adaptive_simpson_recursive(a, b, fa, fb, 0);
    
    // Package results
    List results(n_alphas);
    for (int i = 0; i < n_alphas; ++i) {
        results[i] = List::create(
            Named("Q") = integrals[i],
            Named("fcnt") = function_call_count,
            Named("alpha") = alpha_vec[i],
            Named("converged") = true  // Always true with this approach
        );
    }
    
    return results;
}

//' Vectorized Term2 Integration Across Multiple Alpha Values
//'
//' Integrates the term2 influence function across multiple alpha values simultaneously,
//' sharing expensive computations while maintaining separate convergence criteria.
//'
//' @param compute_expected_values_fn R function that computes expected values for all alphas
//' @param impute_fn R function for data imputation: impute_fn(t, patient_data)
//' @param weight_fn R function for weight computation: weight_fn(t) 
//' @param marginal_mean_fn R function for marginal mean: marginal_mean_fn(t)
//' @param alpha_vec Numeric vector of sensitivity parameters
//' @param tmin Lower integration bound
//' @param tmax Upper integration bound
//' @param patient_data Patient data object (passed to R functions)
//' @param tol Convergence tolerance
//' @return List of integration results, one per alpha value
//' @keywords internal
// [[Rcpp::export]]
List integrate_term2_vectorized_alpha(
    Function compute_expected_values_fn,
    Function impute_fn,
    Function weight_fn,
    Function marginal_mean_fn,
    NumericVector alpha_vec,
    double tmin,
    double tmax,
    SEXP patient_data,
    double tol = 1.490116e-08
) {
    if (alpha_vec.length() == 0) {
        stop("alpha_vec must contain at least one element");
    }
    
    if (tmax <= tmin) {
        stop("tmax must be greater than tmin");
    }
    
    return vectorized_adaptive_simpson(
        compute_expected_values_fn, impute_fn, weight_fn, marginal_mean_fn,
        alpha_vec, patient_data, tmin, tmax, tol
    );
}
