#' Simulate SensIAT Data
#'
#' Generates simulated longitudinal data following the SensIAT model structure,
#' alternating between observation time generation (Cox proportional hazards model)
#' and outcome generation.
#'
#' @param n_subjects Number of subjects to simulate.
#' @param End Maximum follow-up time.
#' @param intensity_coef Coefficient for the effect of previous outcome on observation intensity.
#'        Can be a scalar or vector (one per visit number stratum).
#' @param outcome_coef Named list of coefficients for outcome model including:
#'   * `prev_outcome` - coefficients for natural spline of previous outcome (length 3)
#'   * `time` - coefficient for time since baseline
#'   * `delta_time` - coefficient for time since last observation
#'   * `intercept` - intercept term
#' @param baseline_hazard Baseline hazard function. Either a function of time and visit number,
#'        or a numeric value for constant baseline hazard.
#' @param outcome_sd Standard deviation of the outcome residuals.
#' @param baseline_outcome_fn Optional function to generate baseline outcome value for each subject.
#' @param initial_outcome_mean Mean of the initial (baseline) outcome.
#' @param initial_outcome_sd Standard deviation of the initial outcome.
#' @param max_visits Maximum number of visits per subject (to prevent infinite loops).
#' @param seed Random seed for reproducibility.
#' @param link Link function for outcome model. One of "identity", "log", or "logit".
#'        Determines the scale on which the outcome model operates.
#' @param intensity_fn Optional function to compute intensity (hazard) of observation.
#'        If provided, should take arguments (time, prev_outcome, visit_num) and return
#'        a scalar intensity value. If NULL (default), intensity is computed from
#'        intensity_coef and baseline_hazard.
#' @param intensity_bound Upper bound on intensity for rejection sampling.
#'        Required if intensity_fn is provided. Represents the supremum of the intensity
#'        function on the interval of interest.
#' @param outcome_model Optional fitted single-index outcome model. If provided,
#'        outcomes for follow-up visits are generated from the fitted model via
#'        [make_single_index_simulator()].
#' @param outcome_simulator Optional simulator function for follow-up outcomes.
#'        When provided, it overrides the internal outcome generation function.
#'        This function should accept `prev_outcome`, `time`, `delta_time`, and
#'        optionally `newdata`.
#'
#' @return A tibble with columns:
#'   * `Subject_ID` - Subject identifier
#'   * `Time` - Observation time
#'   * `Outcome` - Observed outcome value
#'   * Additional columns may be added for internal use
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Default usage (uses exponential gaps derived from intensity_coef and baseline_hazard)
#' sim_data <- simulate_SensIAT_data(
#'     n_subjects = 100,
#'     End = 830,
#'     intensity_coef = -0.05,
#'     outcome_coef = list(
#'         prev_outcome = c(0.7, -0.1, 0.05),
#'         time = -0.001,
#'         delta_time = -0.002,
#'         intercept = 2
#'     ),
#'     baseline_hazard = 0.005,
#'     outcome_sd = 1.5
#' )
#'
#' # Example with custom intensity function and thinning (Exp(lambda_star))
#' intensity_fn <- function(t, prev_outcome, visit_num) {
#'   lambda0 <- 0.005
#'   gamma  <- -0.05
#'   lambda0 * exp(gamma * prev_outcome)
#' }
#' sim_data2 <- simulate_SensIAT_data(
#'     n_subjects = 50,
#'     End = 200,
#'     seed = 123,
#'     intensity_fn = intensity_fn,
#'     intensity_bound = 0.05,
#'     max_visits = 20
#' )
#'
#' # Example using a fitted single-index outcome model to generate outcomes
#' # via the fitted conditional distribution.
#' #
#' # Note: this example uses a fitted model object and is intended for
#' # parametric bootstrap-style simulation.
#' #
#' # 
#' # 
#' # outcome_model <- fit_SensIAT_single_index_fixed_coef_model(
#' #     Outcome ~ prev_outcome + time + delta_time,
#' #     data = training_data,
#' #     id = Subject_ID
#' # )
#' # sim_data3 <- simulate_SensIAT_data(
#' #     n_subjects = 50,
#' #     End = 200,
#' #     seed = 123,
#' #     outcome_model = outcome_model,
#' #     intensity_fn = intensity_fn,
#' #     intensity_bound = 0.05,
#' #     max_visits = 20
#' # )
#' }
simulate_SensIAT_data <- function(n_subjects,
                                  End,
                                  intensity_coef = -0.05,
                                  outcome_coef = list(
                                      prev_outcome = c(0.7, -0.1, 0.05),
                                      time = -0.001,
                                      delta_time = -0.002,
                                      intercept = 2
                                  ),
                                  baseline_hazard = 0.005,
                                  outcome_sd = 1.5,
                                  baseline_outcome_fn = NULL,
                                  initial_outcome_mean = 5,
                                  initial_outcome_sd = 2,
                                  max_visits = 50,
                                  seed = NULL,
                                  link = "identity",
                                  intensity_fn = NULL,
                                  intensity_bound = NULL,
                                  outcome_model = NULL,
                                  outcome_simulator = NULL) {
    # Input validation
    if (!is.numeric(n_subjects) || n_subjects < 1 || n_subjects != floor(n_subjects)) {
        stop("n_subjects must be a positive integer")
    }
    if (!is.numeric(End) || End <= 0) {
        stop("End must be a positive number")
    }
    if (!is.numeric(max_visits) || max_visits < 2) {
        stop("max_visits must be at least 2")
    }
    link <- match.arg(link, choices = c("identity", "log", "logit"))
    
    # Validate intensity function parameters
    if (!is.null(intensity_fn)) {
        if (!is.function(intensity_fn)) {
            stop("intensity_fn must be a function")
        }
        if (is.null(intensity_bound)) {
            stop("intensity_bound must be provided when using intensity_fn")
        }
        if (!is.numeric(intensity_bound) || intensity_bound <= 0) {
            stop("intensity_bound must be a positive number")
        }
    }

    # Validate optional outcome simulator inputs
    if (!is.null(outcome_model) && !is.null(outcome_simulator)) {
        stop("Cannot provide both outcome_model and outcome_simulator")
    }
    if (!is.null(outcome_model)) {
        outcome_simulator <- make_single_index_simulator(outcome_model)
    }
    if (!is.null(outcome_simulator) && !is.function(outcome_simulator)) {
        stop("outcome_simulator must be a function")
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # Simulate data for each subject
    results <- purrr::map_dfr(
        seq_len(n_subjects),
        ~ {
            # Keep trying until we get a subject with at least one follow-up
            repeat {
                subject_data <- simulate_single_subject(
                    subject_id = .x,
                    End = End,
                    intensity_coef = intensity_coef,
                    outcome_coef = outcome_coef,
                    baseline_hazard = baseline_hazard,
                    outcome_sd = outcome_sd,
                    baseline_outcome_fn = baseline_outcome_fn,
                    initial_outcome_mean = initial_outcome_mean,
                    initial_outcome_sd = initial_outcome_sd,
                    max_visits = max_visits,
                    link = link,
                    intensity_fn = intensity_fn,
                    intensity_bound = intensity_bound,
                    outcome_simulator = outcome_simulator
                )
                # Check if subject has at least one follow-up observation
                if (max(subject_data$Time) > 0) {
                    break
                }
            }
            subject_data
        }
    )
    
    results
}

#' Simulate Data for a Single Subject
#'
#' Internal function to simulate observation times and outcomes for one subject,
#' alternating between the two processes.
#'
#' @inheritParams simulate_SensIAT_data
#' @param subject_id ID for this subject.
#'
#' @return A tibble with the subject's longitudinal data.
#' @keywords internal
simulate_single_subject <- function(subject_id,
                                    End,
                                    intensity_coef,
                                    outcome_coef,
                                    baseline_hazard,
                                    outcome_sd,
                                    max_visits,
                                    baseline_outcome_fn = NULL,
                                    initial_outcome_mean,
                                    initial_outcome_sd,
                                    link = "identity",
                                    intensity_fn = NULL,
                                    intensity_bound = NULL,
                                    outcome_simulator = NULL) {
    # Initialize vectors to store visit data
    times <- numeric(max_visits)
    outcomes <- numeric(max_visits)
    
    # Baseline observation (time = 0)
    times[1] <- 0
    if (!is.null(baseline_outcome_fn)) {
        if (is.function(baseline_outcome_fn)) {
            baseline_outcome <- baseline_outcome_fn()
            if (!is.numeric(baseline_outcome) || length(baseline_outcome) != 1 || is.na(baseline_outcome)) {
                stop("baseline_outcome_fn must return a single numeric value")
            }
        } else if (is.numeric(baseline_outcome_fn)) {
            baseline_outcome <- sample(baseline_outcome_fn, size=1)
        } else {
            stop("baseline_outcome_fn must be a function or a numeric vector.")
        }
    } else {
        if (link == "logit") {
            # For logit link, generate probability then binary outcome
            outcomes[1] <- stats::rbinom(1, size = 1, prob = initial_outcome_mean)
        } else if (link == "log") {
            # For log link, generate count outcome
            outcomes[1] <- stats::rpois(1, lambda = initial_outcome_mean)
        } else {
            # For identity link, generate continuous outcome
            outcomes[1] <- stats::rnorm(1, mean = initial_outcome_mean, sd = initial_outcome_sd)
        }
        outcomes[1] <- stats::rnorm(1, mean = initial_outcome_mean, sd = initial_outcome_sd)
    }
    
    visit_num <- 1
    current_time <- 0
    current_outcome <- outcomes[1]
    
    # Generate subsequent observations
    while (visit_num < max_visits) {
        # Generate next observation time
        next_time <- generate_next_observation_time(
            current_time = current_time,
            current_outcome = current_outcome,
            visit_num = visit_num,
            intensity_coef = intensity_coef,
            baseline_hazard = baseline_hazard,
            End = End,
            intensity_fn = intensity_fn,
            intensity_bound = intensity_bound
        )
        
        # Stop if next observation exceeds End time
        if (next_time >= End) {
            break
        }
        
        # Generate outcome at this observation
        visit_num <- visit_num + 1
        delta_time <- next_time - current_time
        
        if (!is.null(outcome_simulator)) {
            next_outcome <- outcome_simulator(
                prev_outcome = current_outcome,
                time = next_time,
                delta_time = delta_time
            )
            if (!is.numeric(next_outcome) || length(next_outcome) != 1 || is.na(next_outcome)) {
                stop("outcome_simulator must return a single numeric value")
            }
        } else {
            next_outcome <- generate_outcome(
                prev_outcome = current_outcome,
                time = next_time,
                delta_time = delta_time,
                outcome_coef = outcome_coef,
                outcome_sd = outcome_sd,
                link = link
            )
        }
        
        # Store generated visit
        times[visit_num] <- next_time
        outcomes[visit_num] <- next_outcome

        # Update current state
        current_time <- next_time
        current_outcome <- next_outcome
    }

    # Trim to actual number of visits
    n_visits <- visit_num

    # Return as tibble
    tibble::tibble(
        Subject_ID = subject_id,
        Time = times[1:n_visits],
        Outcome = outcomes[1:n_visits]
    )
}

#' Generate Next Observation Time
#'
#' Generate the time of the next observation using a Cox proportional hazards model.
#' The hazard depends on the previous outcome value and is stratified by visit number.
#'
#' @inheritParams simulate_single_subject
#' @param current_time Current time point.
#' @param current_outcome Current outcome value.
#' @param visit_num Current visit number (for stratification).
#'
#' @return Time of next observation.
#' @keywords internal
generate_next_observation_time <- function(current_time,
                                           current_outcome,
                                           visit_num,
                                           intensity_coef = NULL,
                                           baseline_hazard = NULL,
                                           End,
                                           intensity_fn = NULL,
                                           intensity_bound = NULL) {
    # Use thinning (Exp(lambda_star)) rejection sampling if intensity_fn is provided
    if (!is.null(intensity_fn)) {
        # intensity_bound is interpreted as lambda^* in the algorithm
        lambda_star <- intensity_bound
        if (!is.numeric(lambda_star) || lambda_star <= 0) {
            stop("intensity_bound (lambda_star) must be a positive number")
        }

        repeat {
            # Draw candidate gap time from Exp(lambda_star)
            s_candidate <- stats::rexp(1, rate = lambda_star)
            t_candidate <- current_time + s_candidate

            # If candidate time exceeds End, return it (calling code will stop)
            if (t_candidate > End) {
                return(t_candidate)
            }

            # Evaluate intensity at candidate time
            lambda_t <- intensity_fn(t_candidate, current_outcome, visit_num)

            if (!is.numeric(lambda_t) || length(lambda_t) != 1) {
                stop("intensity_fn must return a single numeric value")
            }
            if (lambda_t < 0) {
                stop("Intensity function returned negative value at t=", t_candidate)
            }
            if (lambda_t > lambda_star) {
                stop(sprintf("Intensity function exceeded provided lambda_star at t=%s: intensity=%s, lambda_star=%s",
                             format(t_candidate, digits = 6), format(lambda_t, digits = 6), format(lambda_star, digits = 6)))
            }

            # Accept with probability lambda(t)/lambda_star
            if (stats::runif(1) < lambda_t / lambda_star) {
                return(t_candidate)
            }

            # If rejected, set start time to t_candidate and repeat
            current_time <- t_candidate
            # continue loop
        }
    }
    
    # Fallback: compute intensity from coefficients (backward compatibility)
    # Compute hazard ratio based on current outcome
    # Using non-linear transformation (could use splines in more complex version)
    covariate_effect <- current_outcome
    
    # Get coefficient for this visit stratum
    if (length(intensity_coef) > 1) {
        coef <- intensity_coef[min(visit_num, length(intensity_coef))]
    } else {
        coef <- intensity_coef
    }
    
    hazard_ratio <- exp(coef * covariate_effect)
    
    # Get baseline hazard
    if (is.function(baseline_hazard)) {
        lambda0 <- baseline_hazard(current_time, visit_num)
    } else {
        lambda0 <- baseline_hazard
    }
    
    # Effective hazard
    lambda <- lambda0 * hazard_ratio
    
    # Generate time to next event from exponential distribution
    time_to_event <- stats::rexp(1, rate = lambda)
    
    # Next observation time
    next_time <- current_time + time_to_event
    
    return(next_time)
}

#' Generate Outcome Value
#'
#' Generate an outcome value based on previous outcome, time since baseline,
#' and time since last observation using a non-linear model.
#'
#' The non-linear transformation of previous outcome approximates the natural spline
#' basis (df=3) used in fitting functions. While fitting uses `splines::ns()`,
#' simulation uses a polynomial approximation for computational efficiency during
#' one-at-a-time data generation.
#'
#' @inheritParams simulate_single_subject
#' @param prev_outcome Previous outcome value.
#' @param time Current time (time since baseline).
#' @param delta_time Time since last observation.
#'
#' @return Generated outcome value.
#' @keywords internal
generate_outcome <- function(prev_outcome,
                             time,
                             delta_time,
                             outcome_coef,
                             outcome_sd,
                             link = "identity") {
    # Non-linear transformation of previous outcome
    # Polynomial approximation of natural spline basis (df=3)
    # This provides similar non-linear behavior to ns(prev_outcome, df=3)
    # used in fitting, while being suitable for single-value generation
    basis1 <- prev_outcome
    basis2 <- prev_outcome^2 / 10
    basis3 <- prev_outcome^3 / 100
    
    # Linear predictor
    eta <- outcome_coef$intercept +
        outcome_coef$prev_outcome[1] * basis1 +
        outcome_coef$prev_outcome[2] * basis2 +
        outcome_coef$prev_outcome[3] * basis3 +
        outcome_coef$time * time +
        outcome_coef$delta_time * delta_time
    
    # Apply link function and generate outcome
    if (link == "logit") {
        # Logit link: generate binary outcome
        prob <- plogis(eta)
        outcome <- stats::rbinom(1, size = 1, prob = prob)
    } else if (link == "log") {
        # Log link: generate count outcome
        lambda <- exp(eta)
        outcome <- stats::rpois(1, lambda = lambda)
    } else {
        # Identity link: generate continuous outcome
        outcome <- stats::rnorm(1, mean = eta, sd = outcome_sd)
    }
    
    return(outcome)
}

#' Simulate Treatment and Control Groups
#'
#' Generate simulated data for both treatment and control groups with
#' potentially different parameters.
#'
#' @inheritParams simulate_SensIAT_data
#' @param treatment_effect Additive treatment effect on outcomes (added to intercept on link scale).
#' @param treatment_intensity_effect Multiplicative effect on observation intensity
#'        (values < 1 mean fewer observations in treatment group).
#' @param link Link function for outcome model. One of "identity", "log", or "logit".
#' @param intensity_fn Optional function to compute intensity (hazard) of observation.
#'        If provided, should take arguments (time, prev_outcome, visit_num) and return
#'        a scalar intensity value. If NULL (default), intensity is computed from
#'        `intensity_coef` and `baseline_hazard`.
#' @param intensity_bound Upper bound on intensity for rejection sampling.
#'        Required if `intensity_fn` is provided. Represents the supremum of the intensity
#'        function on the interval of interest.
#'
#' @return A tibble with an additional `Treatment` column indicating group assignment.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Default treatment/control simulation (uses exponential gaps derived from coefficients)
#' sim_data <- simulate_SensIAT_two_groups(
#'     n_subjects = 100,
#'     End = 830,
#'     treatment_effect = 1.5,
#'     treatment_intensity_effect = 0.9
#' )
#'
#' # Example using custom intensity with thinning
#' intensity_fn <- function(t, prev_outcome, visit_num) {
#'   lambda0 <- 0.005
#'   gamma  <- -0.05
#'   lambda0 * exp(gamma * prev_outcome)
#' }
#' sim_data2 <- simulate_SensIAT_two_groups(
#'     n_subjects = 100,
#'     End = 200,
#'     seed = 123,
#'     intensity_fn = intensity_fn,
#'     intensity_bound = 0.05,
#'     treatment_effect = 1.0,
#'     treatment_intensity_effect = 0.9
#' )
#' }
simulate_SensIAT_two_groups <- function(n_subjects,
                                        End,
                                        intensity_coef = -0.05,
                                        outcome_coef = list(
                                            prev_outcome = c(0.7, -0.1, 0.05),
                                            time = -0.001,
                                            delta_time = -0.002,
                                            intercept = 2
                                        ),
                                        baseline_hazard = 0.005,
                                        outcome_sd = 1.5,
                                        initial_outcome_mean = 5,
                                        initial_outcome_sd = 2,
                                        max_visits = 50,
                                        treatment_effect = 0,
                                        treatment_intensity_effect = 1,
                                        seed = NULL,
                                        link = "identity",
                                        intensity_fn = NULL,
                                        intensity_bound = NULL,
                                        outcome_model = NULL,
                                        outcome_simulator = NULL) {
    # Input validation
    if (!is.numeric(n_subjects) || n_subjects < 2 || n_subjects != floor(n_subjects)) {
        stop("n_subjects must be a positive integer >= 2")
    }
    if (!is.numeric(End) || End <= 0) {
        stop("End must be a positive number")
    }
    link <- match.arg(link, choices = c("identity", "log", "logit"))
    
    # Validate intensity function parameters
    if (!is.null(intensity_fn)) {
        if (!is.function(intensity_fn)) {
            stop("intensity_fn must be a function")
        }
        if (is.null(intensity_bound)) {
            stop("intensity_bound must be provided when using intensity_fn")
        }
        if (!is.numeric(intensity_bound) || intensity_bound <= 0) {
            stop("intensity_bound must be a positive number")
        }
    }

    # Validate optional outcome simulator inputs
    if (!is.null(outcome_model) && !is.null(outcome_simulator)) {
        stop("Cannot provide both outcome_model and outcome_simulator")
    }
    if (!is.null(outcome_model)) {
        outcome_simulator <- make_single_index_simulator(outcome_model)
    }
    if (!is.null(outcome_simulator) && !is.function(outcome_simulator)) {
        stop("outcome_simulator must be a function")
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # Number of subjects per group
    n_control <- floor(n_subjects / 2)
    n_treatment <- n_subjects - n_control
    
    # Simulate control group
    control_data <- simulate_SensIAT_data(
        n_subjects = n_control,
        End = End,
        intensity_coef = intensity_coef,
        outcome_coef = outcome_coef,
        baseline_hazard = baseline_hazard,
        outcome_sd = outcome_sd,
        initial_outcome_mean = initial_outcome_mean,
        initial_outcome_sd = initial_outcome_sd,
        max_visits = max_visits,
        seed = NULL,  # Already set seed above
        link = link,
        intensity_fn = intensity_fn,
        intensity_bound = intensity_bound,
        outcome_model = outcome_model,
        outcome_simulator = outcome_simulator
    ) |>
        dplyr::mutate(Group = "Control")
    
    # Modify parameters for treatment group
    treatment_outcome_coef <- outcome_coef
    treatment_outcome_coef$intercept <- treatment_outcome_coef$intercept + treatment_effect
    
    treatment_baseline_hazard <- baseline_hazard * treatment_intensity_effect
    
    # Simulate treatment group
    treatment_data <- simulate_SensIAT_data(
        n_subjects = n_treatment,
        End = End,
        intensity_coef = intensity_coef,
        outcome_coef = treatment_outcome_coef,
        baseline_hazard = treatment_baseline_hazard,
        outcome_sd = outcome_sd,
        initial_outcome_mean = initial_outcome_mean,
        initial_outcome_sd = initial_outcome_sd,
        max_visits = max_visits,
        seed = NULL,
        link = link,
        intensity_fn = intensity_fn,
        intensity_bound = intensity_bound,
        outcome_model = outcome_model,
        outcome_simulator = outcome_simulator
    ) |>
        dplyr::mutate(
            Subject_ID = .data$Subject_ID + n_control,  # Offset IDs
            Group = "Treatment"
        )
    
    # Combine groups
    result <- dplyr::bind_rows(control_data, treatment_data) |>
        dplyr::arrange(.data$Group, .data$Subject_ID, .data$Time)
    
    return(result)
}
