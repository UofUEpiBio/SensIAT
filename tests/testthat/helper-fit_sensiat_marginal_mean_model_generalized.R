# Helper functions for fit_SensIAT_marginal_mean_model_generalized tests

# Generate test data using in-package simulation functions
generate_test_data <- function(link = "identity", n_subjects = 20, seed = 123) {
    # Set outcome parameters based on link function
    if (link == "identity") {
        initial_outcome_mean <- 50
        initial_outcome_sd <- 5
        outcome_coef <- list(
            prev_outcome = c(0.6, -0.1, 0.05),
            time = -0.002,
            delta_time = -0.001,
            intercept = 10
        )
        outcome_sd <- 3
        End <- 150
        knots <- c(30, 90, 150)
    } else if (link == "log") {
        initial_outcome_mean <- 8  # Poisson lambda (higher for more stability)
        initial_outcome_sd <- NULL  # Not used for Poisson
        outcome_coef <- list(
            prev_outcome = c(0.08, 0.02, -0.01),
            time = 0.0005,
            delta_time = 0.001,
            intercept = 1.8
        )
        outcome_sd <- NULL  # Not used for Poisson
        End <- 120
        knots <- c(30, 60, 90)
    } else if (link == "logit") {
        initial_outcome_mean <- 0.5  # Probability for binary (closer to 0.5 for stability)
        initial_outcome_sd <- NULL  # Not used for binary
        outcome_coef <- list(
            prev_outcome = c(0.2, 0.05, -0.02),  # Smaller coefficients
            time = 0.001,
            delta_time = 0.002,
            intercept = 0  # Centered at 0 for balanced probabilities
        )
        outcome_sd <- NULL  # Not used for binary
        End <- 120
        knots <- c(30, 60, 90)
    } else {
        stop("Invalid link function")
    }
    
    # Generate simulated data
    data <- simulate_SensIAT_data(
        n_subjects = n_subjects,
        End = End,
        intensity_coef = -0.05,
        outcome_coef = outcome_coef,
        baseline_hazard = 0.01,
        outcome_sd = outcome_sd,
        initial_outcome_mean = initial_outcome_mean,
        initial_outcome_sd = initial_outcome_sd,
        max_visits = 20,
        seed = seed,
        link = link
    )
    
    # Add terminal observations
    data_prepared <- prepare_SensIAT_data(
        data,
        id.var = Subject_ID,
        time.var = Time,
        outcome.var = Outcome,
        End = End,
        add.terminal.observations = TRUE
    )
    
    # Fit intensity model
    intensity.model <- survival::coxph(
        survival::Surv(..prev_time.., ..time.., !is.na(..outcome..)) ~ ..prev_outcome..,
        data = data_prepared |> dplyr::filter(..time.. > 0, ..time.. > ..prev_time..)
    )
    
    # Add required attributes for intensity model
    attr(intensity.model, "bandwidth") <- NULL
    attr(intensity.model, "kernel") <- \(x) 0.75 * (1 - (x)**2) * (abs(x) < 1)
    
    # Fit outcome model
    outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
        ..outcome.. ~ ..prev_outcome.. + ..delta_time.. - 1,
        id = ..id..,
        data = data_prepared |> dplyr::filter(..time.. > 0, !is.na(..outcome..))
    )
    
    list(
        data = data_prepared,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        knots = knots,
        End = End
    )
}

# Helper function to create impute_data function
create_impute_fn <- function() {
    function(t, df) {
        data_wl <- df |>
            dplyr::mutate(
                ..prev_time.. = .data$..time..,
                ..prev_outcome.. = .data$..outcome..,
                ..delta_time.. = 0
            )
        extrapolate_from_last_observation(t, data_wl, "..time..", slopes = c("..delta_time.." = 1))
    }
}
