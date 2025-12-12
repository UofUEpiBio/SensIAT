#include "vectorized_integration.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

// Structure to track integration state for each alpha
struct AlphaIntegrationState {
    int alpha_index;
    double current_estimate;
    double current_error;
    bool converged;
    int function_calls;
    
    AlphaIntegrationState(int idx) : 
        alpha_index(idx), current_estimate(0.0), current_error(0.0), 
        converged(false), function_calls(0) {}
};

// Structure to manage overall integration state
struct VectorizedIntegrationState {
    std::vector<AlphaIntegrationState> alpha_states;
    int total_function_calls;
    double tolerance;
    
    VectorizedIntegrationState(int n_alphas, double tol) : 
        total_function_calls(0), tolerance(tol) {
        alpha_states.reserve(n_alphas);
        for (int i = 0; i < n_alphas; ++i) {
            alpha_states.emplace_back(i);
        }
    }
    
    bool all_converged() const {
        return std::all_of(alpha_states.begin(), alpha_states.end(),
                          [](const AlphaIntegrationState& state) { 
                              return state.converged; 
                          });
    }
};

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

// Adaptive Simpson's rule for vectorized integration
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
    VectorizedIntegrationState integration_state(n_alphas, tolerance);
    
    // Initial evaluation points
    double h = 0.13579 * (b - a);
    std::vector<double> x_vals = {a, a + h, a + 2*h, 0.5*(a + b), b - 2*h, b - h, b};
    std::vector<NumericMatrix> y_vals;
    
    // Evaluate integrand at initial points
    for (double x : x_vals) {
        y_vals.push_back(evaluate_integrand_all_alphas(
            x, compute_expected_values_fn, impute_fn, weight_fn, 
            marginal_mean_fn, alpha_vec, patient_data
        ));
        integration_state.total_function_calls++;
    }
    
    int weight_length = y_vals[0].nrow();
    std::vector<NumericVector> final_integrals(n_alphas);
    for (int i = 0; i < n_alphas; ++i) {
        final_integrals[i] = NumericVector(weight_length, 0.0);
    }
    
    // Recursive adaptive integration helper
    std::function<void(int, int, int, int, int, int)> adaptive_helper;
    adaptive_helper = [&](int idx_a, int idx_c, int idx_e, int eval_b, int eval_d, int depth) {
        
        if (integration_state.total_function_calls > 10000) {
            stop("Too many function evaluations; stopping integration");
        }
        
        if (depth > 20) {
            // Force convergence at maximum depth
            for (auto& state : integration_state.alpha_states) {
                state.converged = true;
            }
            return;
        }
        
        double xa = x_vals[idx_a], xc = x_vals[idx_c], xe = x_vals[idx_e];
        double h_seg = xe - xa;
        
        if (h_seg < 1e-12) {
            return; // Segment too small
        }
        
        double xb = 0.5 * (xa + xc);
        double xd = 0.5 * (xc + xe);
        
        // Evaluate at midpoints if needed
        NumericMatrix fb, fd;
        if (eval_b < 0) {
            fb = evaluate_integrand_all_alphas(
                xb, compute_expected_values_fn, impute_fn, weight_fn,
                marginal_mean_fn, alpha_vec, patient_data
            );
            x_vals.push_back(xb);
            y_vals.push_back(fb);
            eval_b = x_vals.size() - 1;
            integration_state.total_function_calls++;
        } else {
            fb = y_vals[eval_b];
        }
        
        if (eval_d < 0) {
            fd = evaluate_integrand_all_alphas(
                xd, compute_expected_values_fn, impute_fn, weight_fn,
                marginal_mean_fn, alpha_vec, patient_data
            );
            x_vals.push_back(xd);
            y_vals.push_back(fd);
            eval_d = x_vals.size() - 1;
            integration_state.total_function_calls++;
        } else {
            fd = y_vals[eval_d];
        }
        
        // Simpson's rule estimates for each alpha
        std::vector<bool> alpha_converged(n_alphas, false);
        
        for (int alpha_idx = 0; alpha_idx < n_alphas; ++alpha_idx) {
            if (integration_state.alpha_states[alpha_idx].converged) {
                alpha_converged[alpha_idx] = true;
                continue;
            }
            
            // Simpson's 1/3 rule (3 points): (h/6)(f(a) + 4f(c) + f(e))
            NumericVector Q1(weight_length);
            for (int j = 0; j < weight_length; ++j) {
                Q1[j] = (h_seg / 6.0) * (
                    y_vals[idx_a](j, alpha_idx) + 
                    4.0 * y_vals[idx_c](j, alpha_idx) + 
                    y_vals[idx_e](j, alpha_idx)
                );
            }
            
            // Simpson's 1/3 rule (5 points): (h/12)(f(a) + 4f(b) + 2f(c) + 4f(d) + f(e))
            NumericVector Q2(weight_length);
            for (int j = 0; j < weight_length; ++j) {
                Q2[j] = (h_seg / 12.0) * (
                    y_vals[idx_a](j, alpha_idx) + 
                    4.0 * fb(j, alpha_idx) +
                    2.0 * y_vals[idx_c](j, alpha_idx) +
                    4.0 * fd(j, alpha_idx) +
                    y_vals[idx_e](j, alpha_idx)
                );
            }
            
            // Romberg extrapolation: Q = Q2 + (Q2 - Q1)/15
            NumericVector Q(weight_length);
            double max_diff = 0.0;
            for (int j = 0; j < weight_length; ++j) {
                Q[j] = Q2[j] + (Q2[j] - Q1[j]) / 15.0;
                max_diff = std::max(max_diff, std::abs(Q2[j] - Q[j]));
            }
            
            // Check convergence for this alpha
            if (max_diff < tolerance) {
                // Add to final integral and mark as converged
                for (int j = 0; j < weight_length; ++j) {
                    final_integrals[alpha_idx][j] += Q[j];
                }
                integration_state.alpha_states[alpha_idx].converged = true;
                alpha_converged[alpha_idx] = true;
            }
        }
        
        // If not all alphas converged, recurse on sub-intervals
        bool any_unconverged = std::any_of(alpha_converged.begin(), alpha_converged.end(),
                                          [](bool conv) { return !conv; });
        
        if (any_unconverged && depth < 20) {
            adaptive_helper(idx_a, eval_b, idx_c, -1, -1, depth + 1);
            adaptive_helper(idx_c, eval_d, idx_e, -1, -1, depth + 1);
        } else {
            // Force add remaining estimates for unconverged alphas
            for (int alpha_idx = 0; alpha_idx < n_alphas; ++alpha_idx) {
                if (!alpha_converged[alpha_idx]) {
                    NumericVector Q2(weight_length);
                    for (int j = 0; j < weight_length; ++j) {
                        Q2[j] = (h_seg / 12.0) * (
                            y_vals[idx_a](j, alpha_idx) + 
                            4.0 * fb(j, alpha_idx) +
                            2.0 * y_vals[idx_c](j, alpha_idx) +
                            4.0 * fd(j, alpha_idx) +
                            y_vals[idx_e](j, alpha_idx)
                        );
                    }
                    for (int j = 0; j < weight_length; ++j) {
                        final_integrals[alpha_idx][j] += Q2[j];
                    }
                }
            }
        }
    };
    
    // Start recursive integration on initial intervals
    adaptive_helper(0, 1, 2, -1, -1, 0);  // [x0, x1, x2]
    adaptive_helper(2, 3, 4, -1, -1, 0);  // [x2, x3, x4] 
    adaptive_helper(4, 5, 6, -1, -1, 0);  // [x4, x5, x6]
    
    // Package results
    List results(n_alphas);
    for (int i = 0; i < n_alphas; ++i) {
        results[i] = List::create(
            Named("Q") = final_integrals[i],
            Named("fcnt") = integration_state.total_function_calls,
            Named("alpha") = alpha_vec[i],
            Named("converged") = integration_state.alpha_states[i].converged
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
//' @export
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