#' @describeIn compute_influence_terms Method for Generalized Linear Models (GLM).
#' @export
compute_influence_terms.glm <-
    function(outcome.model,
             intensity.model,
             alpha,
             data,
             base,
             tolerance = .Machine$double.eps^(1 / 3),
             na.action = na.fail,
             id = NULL,
             time = NULL,
             ...) {
        if (!missing(id)) {
            id <- rlang::ensym(id)
        }
        if (missing(id)) {
            # Try to extract from data or use a default
            if ("..id.." %in% names(data)) {
                id <- rlang::sym("..id..")
            } else {
                rlang::abort("id must be provided for GLM influence terms")
            }
        }
        if (!missing(time)) {
            time <- rlang::ensym(time)
        }
        if (missing(time)) {
            time <- rlang::sym(rlang::f_lhs(terms(intensity.model))[[3]])
        }

        # Filter to follow-up only (time > 0)
        followup_only <- data |>
            dplyr::filter(0 < !!time) |>
            dplyr::select(
                !!id, !!time,
                dplyr::all_of(all.vars(terms(intensity.model))),
                dplyr::all_of(all.vars(terms(outcome.model)))
            ) |>
            na.omit()

        # Get model frame and matrices
        mf_followup <- model.frame(terms(outcome.model), data = followup_only)
        X_followup <- model.matrix(outcome.model, data = mf_followup)
        Y_followup <- model.response(mf_followup)
        ids_followup <- dplyr::pull(followup_only, !!id)
        times_followup <- dplyr::pull(followup_only, !!time)
        uids <- sort(unique(ids_followup))

        # Estimate baseline intensity
        intensity <-
            estimate_baseline_intensity(
                intensity.model = intensity.model,
                data = followup_only
            )
        exp_gamma <- predict(intensity.model, newdata = followup_only, type = "risk", reference = "zero")
        intensity_weights <- intensity$baseline_intensity * exp_gamma

        # Compute term 1: weighted deviations
        term1 <- compute_glm_influence_term_1_for_all(
            data = followup_only,
            glm_model = outcome.model,
            times_all = times_followup,
            ids_all = ids_followup,
            Y_all = Y_followup,
            alpha = alpha,
            intensity_weights = intensity_weights,
            base = base
        )
        
        term1.by.id <- purrr::map(uids, \(uid) {
            colSums(term1[ids_followup == uid, , drop = FALSE], na.rm = TRUE)
        }) |> do.call(rbind, args = _)

        # Compute term 2: integration term
        outcome_var <- rlang::f_lhs(terms(outcome.model))
        integration_data <-
            data |>
            dplyr::select(-dplyr::any_of("..prev_outcome..")) |>
            dplyr::mutate(
                ..prev_outcome..  = !!outcome_var,
                ..delta_time..    = 0,
                ..prev_time..     = !!time,
            )

        term2 <- compute_glm_influence_term_2_for_all_patients(
            outcome.model = outcome.model,
            integration_data = integration_data,
            alpha = alpha,
            base = base,
            tol = tolerance,
            id = !!id,
            time = !!time,
            ...
        )
        
        # Combine terms
        dplyr::full_join(
            tibble::tibble(id = uids, term1 = term1.by.id),
            tibble::tibble(id = attr(term2, "id"), term2 = term2),
            by = "id"
        ) |>
            dplyr::mutate(
                term2 = ifelse(is.na(term2), 0, term2),
                term1 = ifelse(is.na(term1), 0, term1),
                alpha = alpha
            )
    }


#' Compute GLM Influence Term 1 for All Observations
#'
#' @keywords internal
compute_glm_influence_term_1_for_all <-
    function(data,
             glm_model,
             times_all,
             ids_all,
             Y_all,
             alpha,
             intensity_weights,
             base) {
        
        lower <- base@knots[base@order]
        upper <- base@knots[length(base@knots) - base@order + 1]

        # Create model frame from the data to ensure correct column names
        mf_data <- model.frame(terms(glm_model), data = data, na.action = na.pass)
        
        # Get expected values for each observation
        expected_vals <- compute_SensIAT_expected_values(
            model = glm_model,
            alpha = alpha,
            new.data = mf_data
        )
        
        # Compute weighted deviations
        deviations <- (Y_all - expected_vals$E_Yexp_alphaY / expected_vals$E_exp_alphaY) / 
                      intensity_weights

        # Evaluate spline basis at each time point and multiply by deviation
        term1 <- matrix(NA_real_, nrow = length(times_all), ncol = dim(base)[2])
        for (i in seq_along(times_all)) {
            if ((times_all[i] <= lower) || (times_all[i] >= upper)) {
                term1[i, ] <- 0
                next
            }
            
            basis_val <- pcoriaccel_evaluate_basis(base, times_all[i])
            term1[i, ] <- basis_val * deviations[i]
        }
        
        return(term1)
    }


#' Compute GLM Influence Term 2 for All Patients
#'
#' @keywords internal
compute_glm_influence_term_2_for_all_patients <-
    function(outcome.model,
             integration_data,
             alpha,
             base,
             id, time,
             tol = .Machine$double.eps^(1/3),
             ...) {
        
        id <- rlang::ensym(id)
        time <- rlang::ensym(time)
        
        ids <- dplyr::pull(integration_data, !!id)
        times <- dplyr::pull(integration_data, !!time)
        uids <- unique(ids)
        
        # Compute slope (derivative with respect to time)
        slope <- compute_slope(
            outcome.model = outcome.model,
            time.vars = c(rlang::as_string(time)),
            ...
        )

        term2 <- matrix(NA_real_, nrow = length(uids), ncol = dim(base)[2])
        
        for (i in seq_along(uids)) {
            patient_indices <- which(ids == uids[i])
            patient_data <- integration_data[patient_indices, , drop = FALSE]
            patient_times <- times[patient_indices]
            
            result <- compute_glm_influence_term_2_for_individual(
                patient_data = patient_data,
                patient_times = patient_times,
                outcome.model = outcome.model,
                alpha = alpha,
                base = base,
                slope = slope,
                tol = tol
            )
            
            term2[i, ] <- result
        }
        
        structure(term2, id = uids)
    }


#' Compute GLM Influence Term 2 for Individual Patient
#'
#' @keywords internal
compute_glm_influence_term_2_for_individual <-
    function(patient_data,
             patient_times,
             outcome.model,
             alpha,
             base,
             slope,
             tol = .Machine$double.eps^(1/3)) {
        
        lower <- base@knots[base@order]
        upper <- base@knots[length(base@knots) - base@order + 1]
        
        # Initialize result
        result <- numeric(dim(base)[2])
        
        # For each observation time for this patient
        for (j in seq_along(patient_times)) {
            t_j <- patient_times[j]
            
            if (t_j <= lower || t_j >= upper) {
                next
            }
            
            # Get expected values at this time point
            obs_data <- patient_data[j, , drop = FALSE]
            expected_vals <- compute_SensIAT_expected_values(
                model = outcome.model,
                alpha = alpha,
                new.data = obs_data
            )
            
            E_Yexp_alphaY <- expected_vals$E_Yexp_alphaY
            E_exp_alphaY <- expected_vals$E_exp_alphaY
            
            # Compute integrand value
            integrand_val <- E_Yexp_alphaY / E_exp_alphaY * slope
            
            # Evaluate basis at this time
            basis_val <- pcoriaccel_evaluate_basis(base, t_j)
            
            # Add contribution (using simple summation as approximation to integral)
            # In practice, this would use proper numerical integration
            result <- result + basis_val * integrand_val
        }
        
        return(result)
    }
