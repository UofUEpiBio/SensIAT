## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 7,
    fig.height = 5
)


## ----setup--------------------------------------------------------------------
library(SensIAT)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(splines)
set.seed(12345)


## ----generate_random_model_for_sim--------------------------------------------
n_subjects <- 25
end_time   <- 365  # 1 year keeps computation tractable

# --- Outcome simulator (log link) -------------------------------------------
# Properties:
#   - Next outcome inversely proportional to previous (beta_prev < 0 on log scale)
#   - Increases slowly with inter-visit gap (sqrt(delta_time) keeps it bounded)
#   - Increases gradually over study time (small beta_time)
#   - Log-normal generation keeps Y > 0 and prevents blow-up
#   - |beta_prev| < 1 guarantees AR(1) stationarity on log scale
#
# Equilibrium: log(Y*) = intercept / (1 - beta_prev)

beta_prev <- -0.35   # inverse proportionality; |b|<1 => stable mean reversion
beta_dt   <-  0.003  # slow sqrt-dampened increase with inter-visit gap
beta_time <-  3e-4   # gradual drift upward over study period
intercept <-  1.2    # log-scale intercept; exp(1.2/(1+0.35)) ≈ 2.4 equilibrium
sigma_eps <-  0.25   # residual SD on log scale

outcome_simulator <- function(prev_outcome, time, delta_time, ...) {
    log_prev <- log(max(prev_outcome, 1e-6))
    eta <- intercept +
        beta_prev * log_prev +
        beta_dt   * sqrt(delta_time) +
        beta_time * time
    exp(stats::rnorm(1, mean = eta, sd = sigma_eps))
}

# Baseline: log-normal so Y_0 > 0; median ≈ 3
baseline_outcome_fn <- function() stats::rlnorm(1, meanlog = log(3), sdlog = 0.4)

# --- Intensity function (Weibull baseline + outcome effect) ------------------
# Why Weibull here?
#   - shape > 1 gives low hazard early and higher hazard later, which helps avoid
#     very short early follow-up gaps.
#   - outcome effect is multiplicative and positive.
#   - hard bounds keep rates in a realistic range.

# Weibull baseline hazard in calendar time t:
#   h0(t) = (shape / scale) * (t / scale)^(shape - 1)
wb_shape <- 1.35
wb_scale <- 220

# Outcome dependence and caps
gamma_y    <- 0.20      # mild positive effect of previous outcome
y_ref      <- 3         # reference outcome
lambda_min <- 1 / 240   # avoid near-zero hazard (very sparse follow-up)
lambda_max <- 1 / 35    # avoid overly short expected gaps

weibull_baseline_hazard <- function(t, shape = wb_shape, scale = wb_scale) {
    tt <- max(t, 1e-6)
    (shape / scale) * (tt / scale)^(shape - 1)
}

intensity_fn <- function(t, prev_outcome, visit_num) {
    h0 <- weibull_baseline_hazard(t)
    outcome_mult <- exp(gamma_y * (log(max(prev_outcome, 1e-6)) - log(y_ref)))
    rate <- h0 * outcome_mult
    min(max(rate, lambda_min), lambda_max)
}
intensity_bound <- lambda_max

cat(sprintf(
    "Outcome: log(Y_next) ~ N(%.2f + (%.2f)*log(Y_prev) + %.3f*sqrt(dt) + %.1e*t,  sd=%.2f)\n",
    intercept, beta_prev, beta_dt, beta_time, sigma_eps
))
cat(sprintf("Equilibrium Y* ≈ %.2f\n", exp(intercept / (1 - beta_prev))))
cat(sprintf(
    "Intensity: Weibull(shape=%.2f, scale=%.0f), bounded to [1/%.0f, 1/%.0f] per day\n",
    wb_shape, wb_scale, 1 / lambda_min, 1 / lambda_max
))

## ----generation_function------------------------------------------------------
# Show E[Next Outcome | Previous Outcome] on log-log axes at three study times
prev_grid  <- exp(seq(log(0.3), log(20), length.out = 300))
times_show <- c(0, 365, 730)

lapply(times_show, function(t) {
    eta <- intercept +
        beta_prev * log(prev_grid) +
        beta_dt   * sqrt(90) +   # typical 90-day gap
        beta_time * t
    data.frame(prev_outcome = prev_grid, expected_next = exp(eta), day = t)
}) %>%
    bind_rows() %>%
    mutate(day = factor(day, labels = paste("Day", times_show))) %>%
    ggplot(aes(x = prev_outcome, y = expected_next, color = day)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
    geom_line(linewidth = 1.1) +
    scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
    scale_color_brewer(palette = "Blues", direction = 1) +
    labs(
        title = "Data Generation: E[Next Outcome | Previous Outcome]",
        subtitle = sprintf(
            "log(Y_next) = %.2f + (%.2f)*log(Y_prev) + %.3f*sqrt(dt) + %.1e*t",
            intercept, beta_prev, beta_dt, beta_time
        ),
        x = "Previous Outcome (log scale)",
        y = "E[Next Outcome] (log scale)",
        color = NULL,
        caption = paste0(
            "Dotted: y = x (no change). Curves below diagonal = mean reversion. ",
            "dt fixed at 90 days. Equilibrium Y* \u2248 ",
            round(exp(intercept / (1 - beta_prev)), 1)
        )
    ) +
    theme_minimal() +
    theme(plot.caption = element_text(hjust = 0))

## ----simulate-----------------------------------------------------------------
sim_data <- simulate_SensIAT_data(
    n_subjects          = n_subjects,
    End                 = end_time,
    baseline_outcome_fn = baseline_outcome_fn,
    outcome_simulator   = outcome_simulator,
    intensity_fn        = intensity_fn,
    intensity_bound     = intensity_bound,
    seed                = 42,
    max_visits          = 30
)

# Display summary of simulated data
cat("Number of subjects:", n_distinct(sim_data$Subject_ID), "\n")
cat("Number of observations:", nrow(sim_data), "\n")
cat("Outcome range:", range(sim_data$Outcome), "\n")
cat("Time range:", range(sim_data$Time), "\n")

head(sim_data, 10)


## ----plot_sim_data------------------------------------------------------------
# Plot outcomes for a random sample of subjects
sample_subjects <- sample(unique(sim_data$Subject_ID), size = 12, replace = FALSE)

sim_data %>%
    filter(Subject_ID %in% sample_subjects) %>%
    ggplot(aes(x = Time, y = Outcome, group = Subject_ID)) +
    geom_line(alpha = 0.5, color = "steelblue") +
    geom_point(alpha = 0.6, color = "steelblue", size = 2) +
    facet_wrap(~Subject_ID, ncol = 4) +
    labs(
        title = "Observed Outcomes over Time",
        subtitle = "Sample of 12 subjects with irregular assessment times",
        x = "Time (days)",
        y = "Outcome",
        caption = "Simulated data with random outcome model"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ----fit_model----------------------------------------------------------------
# Key speed-control parameters for the generalized (log/logit) solver:
#   term2_method  - "fixed_grid" is 2-5x faster than "fast" for multiple alpha values
#   term2_grid_n  - grid points for fixed_grid (fewer = faster, less accurate)
#   alpha         - fewer values = fewer solver iterations
knots <- c(90, 180, 270)  # knots for spline basis in outcome model

model_log <- fit_SensIAT_within_group_model(
    group.data     = sim_data,
    outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
    outcome.args   = list(
        model = ~ ns(..prev_outcome.., df = 3) + scale(..delta_time..) - 1
    ),
    alpha          = 0, #c(-0.3, 0, 0.3),    # start with 3 values while profiling
    id             = "Subject_ID",
    outcome        = "Outcome",
    time           = "Time",
    End            = end_time,
    knots          = knots,
    link           = "log",              # log link -- exercises generalized solver
    loss           = "lp_mse",
    term2_method   = "fixed_grid",       # faster than "fast" for multiple alpha values
    influence.args = list(
        term2_grid_n    = 20,
        BBsolve.control = list(maxit = 200, tol = 1e-4)
    )
)

# Display model summary
print(model_log)


## ----extract_coefs------------------------------------------------------------
coef_matrix <- coef(model_log)
print(coef_matrix)

# Point-estimate marginal mean curve
autoplot(model_log) +
    labs(
        title = "Fitted Marginal Mean (log link)",
        subtitle = sprintf("n=%d subjects, alpha=0", n_subjects)
    )


## ----jackknife----------------------------------------------------------------
# Compute intensity_bound once globally (only 50 grid points + pre-allocated vector = fast)
# Then pass it to parametric_bootstrap_within_group via simulate_args to avoid recomputation per replicate

cat("\n=== Computing shared intensity bound (coarse grid, once globally) ===\n")
t_bound_start <- proc.time()

# Build parametric simulator for bound computation
intensity_sim_for_bound <- SensIAT:::make_parametric_intensity_simulator(
    model_log$models$intensity,
    SensIAT:::get_bootstrap_coeffs(model_log$models$intensity, sample_coefficients = FALSE)
)

# Helper: compute bound with coarse grid (replicates logic from parametric_bootstrap.R)
compute_bound_coarse <- function(intensity_sim, End, initial_outcome_mean, initial_outcome_sd) {
    times_grid <- unique(c(seq(0, End, length.out = 50), 0, End))
    prev_outcomes <- seq(
        initial_outcome_mean - 3 * initial_outcome_sd,
        initial_outcome_mean + 3 * initial_outcome_sd,
        length.out = 5
    )
    prev_outcomes <- prev_outcomes[is.finite(prev_outcomes)]
    
    evals <- numeric(length(times_grid) * length(prev_outcomes) * 5)
    idx <- 0L
    for (tt in times_grid) {
        for (pv in prev_outcomes) {
            for (vn in 1:5) {
                idx <- idx + 1L
                val <- tryCatch(intensity_sim(tt, pv, vn), error = function(e) 0)
                if (!is.finite(val) || val < 0) val <- 0
                evals[idx] <- val
            }
        }
    }
    evals <- evals[seq_len(idx)]
    max(1e-6, max(evals, na.rm = TRUE) * 10000, 1e5)
}

shared_intensity_bound <- compute_bound_coarse(
    intensity_sim_for_bound, 
    End = end_time,
    initial_outcome_mean = mean(sim_data$Outcome[sim_data$Time == 0], na.rm = TRUE),
    initial_outcome_sd = sd(sim_data$Outcome[sim_data$Time == 0], na.rm = TRUE)
)
t_bound_elapsed <- (proc.time() - t_bound_start)[[3]]

cat(sprintf("Shared intensity_bound computed: %.6f (elapsed: %.2f s)\n\n", shared_intensity_bound, t_bound_elapsed))

# Run bootstrap with pre-computed bound (avoids recomputation per replicate)
pbs_results <- parametric_bootstrap_within_group(
    within_group_model  = model_log,
    nboot   = 20,
    seed     = 12345,
    progress = interactive(),
    verbosity = "detailed",
    # NOTE: Outcome simulator is now included (optimized with pre-computed lp0).
    # Run single_index_simulator_benchmark.R to see before/after performance.
    simulate_args = list(
        intensity_bound = shared_intensity_bound,
        outcome_simulator = SensIAT:::make_parametric_single_index_simulator(
            model_log$models$outcome,
            SensIAT:::get_bootstrap_coeffs(model_log$models$outcome, sample_coefficients = FALSE)
        )
    )
)

print(pbs_results)



## ----jackknife_plot------------------------------------------------------------
# autoplot for bootstrap results draws point estimates + 95% CI bands
autoplot(pbs_results) +
    labs(
        title    = "Marginal Mean with 95% Bootstrap CI",
        subtitle = "Log link, parametric bootstrap with coarse grid bound",
        caption  = "Shaded region: bootstrap 95% CI"
    )


