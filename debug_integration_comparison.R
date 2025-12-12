# Debug script to compare integration methods

devtools::load_all()

# Very simple data that should work for both methods
patient_data <- data.frame(
  Time = c(20, 25, 30),
  prev_outcome_lag1 = c(10, 12, 14),
  Outcome = c(12, 14, 16)
)

outcome_model <- lm(Outcome ~ prev_outcome_lag1, data = patient_data)
print("Model summary:")
print(summary(outcome_model))

# Minimal spline setup 
knots <- c(22, 28)
spline.degree <- 1L  # Linear splines
knots_extended <- c(22, knots, 28)

base <- orthogonalsplinebasis::SplineBasis(knots_extended, order = 2L)
n_basis <- dim(base@Matrices)[2]
marginal_beta <- rep(0.05, n_basis)  # Very small coefficients
V_inv <- diag(n_basis)

tmin <- 22
tmax <- 28

# Simple constant impute function
impute_fn <- function(t, data) {
  data.frame(prev_outcome_lag1 = 12)
}

inv_link <- function(x) x

# Test alpha = 0 (should be simplest case)
alpha_val <- 0
print(paste("Testing alpha =", alpha_val))
print(paste("Integration bounds: [", tmin, ",", tmax, "]"))
print(paste("Number of basis functions:", n_basis))

# Original method
cat("Running original method...\n")
orig_result <- compute_term2_influence_original(
  patient_data, outcome_model, base, alpha_val, marginal_beta,
  V_inv, tmin, tmax, impute_fn, inv_link
)
cat("Original result:\n")
print(orig_result)

# Vectorized method
cat("\nRunning vectorized method...\n")
vec_result <- compute_term2_influence_vectorized(
  patient_data, outcome_model, base, c(alpha_val), marginal_beta,
  V_inv, tmin, tmax, impute_fn, inv_link
)
cat("Vectorized result:\n")
print(vec_result$alpha_0$Q)

cat("\nDifference:\n")
diff <- vec_result$alpha_0$Q - orig_result
print(diff)

cat("\nRelative difference:\n")
rel_diff <- abs(diff) / (abs(orig_result) + 1e-8)
print(rel_diff)