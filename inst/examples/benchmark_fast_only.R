# Benchmark FAST integrand only (bypasses Original and Cached)

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(SensIAT)
  }
  library(dplyr)
  library(purrr)
})

# Data and models
data_with_lags <- SensIAT_example_data |>
  dplyr::group_by(Subject_ID) |>
  dplyr::mutate(
    ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
  ) |>
  dplyr::ungroup()

outcome.model <-
  fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
    id = Subject_ID,
    data = data_with_lags |> dplyr::filter(Time > 0)
  )

alpha <- 0
knots <- c(60, 260, 460)
spline.degree <- 3L
knots_extended <- c(
  rep(head(knots, 1), spline.degree),
  knots,
  rep(tail(knots, 1), spline.degree)
)
base <- orthogonalsplinebasis::SplineBasis(knots_extended, order = spline.degree + 1L)
V.inv <- solve(orthogonalsplinebasis::GramMatrix(base))

tmin <- base@knots[base@order]
tmax <- base@knots[length(base@knots) - base@order + 1]

set.seed(123)
observed_outcomes <- data_with_lags$Outcome[!is.na(data_with_lags$Outcome)]
mean_outcome <- mean(observed_outcomes)
beta_test <- rep(log(mean_outcome) / ncol(base), ncol(base)) + rnorm(ncol(base), 0, 0.1)

test_ids <- unique(data_with_lags$Subject_ID)[1:5]
cat("Testing FAST-only with", length(test_ids), "patients\n\n")

make_fast_integrand_for_patient <- function(patient_id) {
  data_i <- data_with_lags[data_with_lags$Subject_ID == patient_id, ]
  f_builder <- tryCatch(
    get("make_term2_integrand_fast", envir = asNamespace("SensIAT")),
    error = function(e) get("make_term2_integrand_fast")
  )
  f_builder(
    outcome.model = outcome.model,
    base = base,
    alpha = alpha,
    patient_times = data_i$Time,
    patient_outcomes = data_i$Outcome,
    marginal_beta = beta_test,
    V_inv = V.inv
  )
}

cat("Running FAST method...\n")
start <- proc.time()[[3]]
results_fast <- map(
  test_ids,
  ~{
    f <- make_fast_integrand_for_patient(.x)
    rslt <- pcoriaccel_integrate_simp(f, tmin, tmax)
    rslt$Q
  }
)
elapsed_fast <- proc.time()[[3]] - start
cat("FAST method completed in", round(elapsed_fast, 2), "seconds\n\n")

orig_stub <- NA_real_
cat("Summary (FAST-only):\n")
cat("  FAST elapsed: ", round(elapsed_fast, 2), "s\n", sep = "")
cat("  Per-patient:  ", round(elapsed_fast / length(test_ids), 3), "s/patient\n", sep = "")
