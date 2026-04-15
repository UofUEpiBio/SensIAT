# Benchmark term2 integration methods
#
# Provides efficient comparison of term2 computation methods by isolating
# just the integration computations from the full model fitting process.

#' Benchmark Term2 Integration Methods
#'
#' Efficiently compares different term2 integration methods by isolating just
#' the term2 computation from the full model fitting. This is useful for
#' performance analysis and method selection without the overhead of repeated
#' model fitting.
#'
#' @param data Data frame with longitudinal observations
#' @param id Unquoted column name for subject identifier
#' @param time Numeric vector of observation times
#' @param outcome.model A fitted outcome model (e.g., from fit_SensIAT_single_index_fixed_coef_model)
#' @param knots Knots for marginal mean model spline basis
#' @param spline_degree Degree of B-spline basis (default: 3)
#' @param alpha Numeric vector of sensitivity parameters to test
#' @param impute_data Function to impute data at arbitrary times: `function(t, patient_data) -> data.frame`
#' @param link Link function: "identity" or "log"
#' @param methods Character vector of methods to benchmark. Options:
#'   "fast", "original", "fixed_grid", "seeded_adaptive", "gauss_legendre"
#' @param grid_sizes For grid-based methods, vector of grid sizes to test (default: c(50, 100, 200))
#' @param n_patients Number of patients to use for benchmark (NULL = all patients)
#' @param n_iterations Number of timing iterations per method
#' @param reference_method Method to use as accuracy reference (default: "fast")
#' @param seed Random seed for reproducibility
#'
#' @return A list with components:
#'   \describe{
#'     \item{timing}{Data frame with timing results per method}
#'     \item{accuracy}{Data frame with accuracy metrics vs reference}
#'     \item{reference_results}{List of reference results per patient}
#'     \item{setup_info}{List with setup parameters and dimensions}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Setup data with lag variables
#' data_with_lags <- SensIAT_example_data |>
#'   dplyr::group_by(Subject_ID) |>
#'   dplyr::mutate(
#'     ..prev_outcome.. = dplyr::lag(Outcome, order_by = Time),
#'     ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
#'     ..delta_time.. = Time - dplyr::lag(Time, order_by = Time)
#'   ) |>
#'   dplyr::ungroup()
#'
#' # Fit outcome model
#' outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
#'   Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
#'   id = Subject_ID,
#'   data = data_with_lags |> dplyr::filter(Time > 0)
#' )
#'
#' # Imputation function
#' impute_fn <- function(t, df) {
#'   extrapolate_from_last_observation(t, df, "Time",
#'     slopes = c("..delta_time.." = 1))
#' }
#'
#' # Run benchmark
#' results <- benchmark_term2_methods(
#'   data = data_with_lags,
#'   id = Subject_ID,
#'   time = data_with_lags$Time,
#'   outcome.model = outcome.model,
#'   knots = c(60, 260, 460),
#'   alpha = c(-0.1, 0, 0.1),
#'   impute_data = impute_fn,
#'   methods = c("fast", "fixed_grid", "seeded_adaptive"),
#'   grid_sizes = c(50, 100),
#'   n_patients = 10
#' )
#'
#' # View timing results
#' print(results$timing)
#' }
benchmark_term2_methods <- function(
    data,
    id,
    time,
    outcome.model,
    knots,
    spline_degree = 3L,
    alpha = 0,
    impute_data,
    link = c("identity", "log"),
    methods = c("fast", "original", "fixed_grid", "seeded_adaptive"),
    grid_sizes = c(50, 100, 200),
    n_patients = NULL,
    n_iterations = 3,
    reference_method = "fast",
    seed = 42
) {
    link <- match.arg(link)
    
    # Validate methods
    valid_methods <- c("fast", "original", "fixed_grid", "seeded_adaptive", "gauss_legendre")
    if (!all(methods %in% valid_methods)) {
        stop("Invalid method(s). Valid options: ", paste(valid_methods, collapse = ", "))
    }
    
    # Get subject IDs
    id_values <- rlang::eval_tidy(rlang::enquo(id), data)
    unique_ids <- unique(id_values)
    
    # Subset patients if requested
    if (!is.null(n_patients) && n_patients < length(unique_ids)) {
        set.seed(seed)
        unique_ids <- sample(unique_ids, n_patients)
    }
    
    # Setup spline basis
    knots_extended <- c(
        rep(head(knots, 1), spline_degree),
        knots,
        rep(tail(knots, 1), spline_degree)
    )
    base <- orthogonalsplinebasis::SplineBasis(knots_extended, order = spline_degree + 1L)
    V <- orthogonalsplinebasis::GramMatrix(base)
    V.inv <- solve(V)
    
    tmin <- base@knots[base@order]
    tmax <- base@knots[length(base@knots) - base@order + 1]
    
    # Initialize beta with reasonable values
    set.seed(seed)
    observed_outcomes <- data$Outcome[!is.na(data$Outcome)]
    mean_outcome <- mean(observed_outcomes)
    if (link == "log") {
        beta_test <- rep(log(mean_outcome) / ncol(base), ncol(base)) + rnorm(ncol(base), 0, 0.1)
    } else {
        beta_test <- rep(mean_outcome / ncol(base), ncol(base)) + rnorm(ncol(base), 0, 0.1)
    }
    
    # Inverse link function
    inv_link <- if (link == "log") exp else identity
    
    # Weight function
    W <- if (link == "log") {
        function(t, beta) {
            B <- as.vector(pcoriaccel_evaluate_basis(base, t))
            mu <- sum(B * beta)
            as.vector((V.inv %*% B) * exp(-mu))
        }
    } else {
        function(t, beta) {
            B <- as.vector(pcoriaccel_evaluate_basis(base, t))
            as.vector(V.inv %*% B)
        }
    }
    
    # Create imputation wrapper that sets up lag variables for term2
    impute_fn_wrapper <- function(t, df) {
        data_wl <- df |>
            dplyr::mutate(
                ..prev_time.. = .data$Time,
                ..prev_outcome.. = .data$Outcome,
                ..delta_time.. = 0
            )
        impute_data(t, data_wl)
    }
    
    # Get patient data subsets
    patient_data_list <- lapply(unique_ids, function(pid) {
        data[id_values == pid, ]
    })
    names(patient_data_list) <- as.character(unique_ids)
    
    # Build method configurations to test
    method_configs <- build_method_configs(methods, grid_sizes)
    
    # Ensure reference method is computed first
    ref_config_idx <- which(sapply(method_configs, function(x) x$label == reference_method))
    if (length(ref_config_idx) == 0) {
        ref_config_idx <- 1
        reference_method <- method_configs[[1]]$label
    }
    
    # ========================================================================
    # BENCHMARK EXECUTION
    # ========================================================================
    
    cat("Benchmarking", length(method_configs), "method configurations\n")
    cat("Patients:", length(unique_ids), "| Alphas:", length(alpha), 
        "| Iterations:", n_iterations, "\n\n")
    
    all_results <- list()
    all_timings <- list()
    
    for (cfg_idx in seq_along(method_configs)) {
        cfg <- method_configs[[cfg_idx]]
        cat("  ", cfg$label, "...", sep = "")
        
        # Collect results and timing for this configuration
        iteration_times <- numeric(n_iterations)
        
        for (iter in seq_len(n_iterations)) {
            t_start <- proc.time()[3]
            
            results_this_iter <- list()
            for (a_idx in seq_along(alpha)) {
                a <- alpha[a_idx]
                patient_results <- list()
                
                for (pid in names(patient_data_list)) {
                    patient_data <- patient_data_list[[pid]]
                    
                    result <- run_term2_method(
                        cfg = cfg,
                        patient_data = patient_data,
                        outcome.model = outcome.model,
                        base = base,
                        alpha = a,
                        marginal_beta = beta_test,
                        V_inv = V.inv,
                        tmin = tmin,
                        tmax = tmax,
                        impute_fn = impute_fn_wrapper,
                        inv_link = inv_link,
                        W = W
                    )
                    patient_results[[pid]] <- result
                }
                results_this_iter[[paste0("alpha_", a)]] <- patient_results
            }
            
            iteration_times[iter] <- proc.time()[3] - t_start
        }
        
        # Store the last iteration's results for accuracy comparison
        all_results[[cfg$label]] <- results_this_iter
        all_timings[[cfg$label]] <- iteration_times
        
        cat(" ", sprintf("%.3fs", mean(iteration_times)), "\n")
    }
    
    # ========================================================================
    # COMPUTE ACCURACY METRICS
    # ========================================================================
    
    cat("\nComputing accuracy metrics vs reference (", reference_method, ")...\n", sep = "")
    
    reference_results <- all_results[[reference_method]]
    accuracy_metrics <- compute_accuracy_metrics(all_results, reference_results, reference_method)
    
    # ========================================================================
    # FORMAT TIMING RESULTS
    # ========================================================================
    
    timing_df <- data.frame(
        method = names(all_timings),
        mean_time = sapply(all_timings, mean),
        sd_time = sapply(all_timings, sd),
        min_time = sapply(all_timings, min),
        max_time = sapply(all_timings, max),
        stringsAsFactors = FALSE
    )
    timing_df$relative_speed <- timing_df$mean_time / min(timing_df$mean_time)
    timing_df <- timing_df[order(timing_df$mean_time), ]
    rownames(timing_df) <- NULL
    
    # Return results
    list(
        timing = timing_df,
        accuracy = accuracy_metrics,
        reference_results = reference_results,
        setup_info = list(
            n_patients = length(unique_ids),
            n_alphas = length(alpha),
            alpha_values = alpha,
            n_basis_functions = ncol(base),
            tmin = tmin,
            tmax = tmax,
            link = link,
            reference_method = reference_method,
            n_iterations = n_iterations
        )
    )
}


# Internal: Build list of method configurations to test
build_method_configs <- function(methods, grid_sizes) {
    configs <- list()
    
    for (method in methods) {
        if (method %in% c("fast", "original")) {
            configs[[length(configs) + 1]] <- list(
                method = method,
                label = method,
                grid_n = NA
            )
        } else if (method %in% c("fixed_grid", "seeded_adaptive", "gauss_legendre")) {
            for (gs in grid_sizes) {
                configs[[length(configs) + 1]] <- list(
                    method = method,
                    label = paste0(method, "_", gs),
                    grid_n = gs
                )
            }
        }
    }
    
    configs
}


# Internal: Run a single term2 method
run_term2_method <- function(cfg, patient_data, outcome.model, base, alpha,
                              marginal_beta, V_inv, tmin, tmax, impute_fn, 
                              inv_link, W) {
    fn <- switch(cfg$method,
        fast = compute_term2_influence_fast,
        original = compute_term2_influence_original,
        fixed_grid = compute_term2_influence_fixed_grid,
        seeded_adaptive = compute_term2_influence_seeded_adaptive,
        gauss_legendre = compute_term2_influence_gauss_legendre
    )
    
    # Build arguments
    args <- list(
        patient_data = patient_data,
        outcome_model = outcome.model,
        base = base,
        alpha = alpha,
        marginal_beta = marginal_beta,
        V_inv = V_inv,
        tmin = tmin,
        tmax = tmax,
        impute_fn = impute_fn,
        inv_link = inv_link,
        W = W
    )
    
    # Add grid_n for grid-based methods
    if (!is.na(cfg$grid_n)) {
        args$n_grid <- cfg$grid_n
    }
    
    do.call(fn, args)
}


# Internal: Compute accuracy metrics comparing results to reference
compute_accuracy_metrics <- function(all_results, reference_results, reference_method) {
    metrics <- list()
    
    for (method_name in names(all_results)) {
        if (method_name == reference_method) next
        
        method_results <- all_results[[method_name]]
        
        # Compute differences across all alphas and patients
        all_diffs <- c()
        
        for (alpha_name in names(reference_results)) {
            ref_alpha <- reference_results[[alpha_name]]
            method_alpha <- method_results[[alpha_name]]
            
            for (pid in names(ref_alpha)) {
                ref_val <- as.vector(ref_alpha[[pid]])
                method_val <- as.vector(method_alpha[[pid]])
                all_diffs <- c(all_diffs, abs(ref_val - method_val))
            }
        }
        
        metrics[[method_name]] <- list(
            max_abs_diff = max(all_diffs),
            mean_abs_diff = mean(all_diffs),
            rmse = sqrt(mean(all_diffs^2)),
            median_abs_diff = median(all_diffs)
        )
    }
    
    # Convert to data frame
    if (length(metrics) == 0) {
        return(data.frame(
            method = character(),
            max_abs_diff = numeric(),
            mean_abs_diff = numeric(),
            rmse = numeric(),
            stringsAsFactors = FALSE
        ))
    }
    
    data.frame(
        method = names(metrics),
        max_abs_diff = sapply(metrics, `[[`, "max_abs_diff"),
        mean_abs_diff = sapply(metrics, `[[`, "mean_abs_diff"),
        rmse = sapply(metrics, `[[`, "rmse"),
        median_abs_diff = sapply(metrics, `[[`, "median_abs_diff"),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}
