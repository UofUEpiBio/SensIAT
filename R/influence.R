


compute_influence_term_1_at_timepoint <- function(
    time,
    outcome,
    prev_outcome,
    Xb_ind, #< individual-level covariates to be passed to pmf_estimator
    pmf_estimator,
    intensity_coef,  #< single coefficient from the intensity model that is the coefficient corresponding to the intensity from the previous outcome.
    alpha, #< sensitivity parameter.
    baseline_lambda, #< baseline intensity for the individual.
    base,
    y = sort(unique(outcome))
){
    if (getOption('PCORI::do_arg_checks', TRUE))
        assert_that(
            rlang::is_scalar_double(time), time > 0,
            rlang::is_scalar_double(outcome),
            rlang::is_scalar_double(prev_outcome),
            rlang::is_scalar_double(Xb_ind),
            is.function(pmf_estimator),
            rlang::is_scalar_double(intensity_coef),
            rlang::is_scalar_double(alpha),
            rlang::is_scalar_double(baseline_lambda),
            is(base, "SplineBasis")
        )

    lower <- base@knots[base@order]
    upper <- base@knots[length(base@knots) - base@order+1]

    if ((time <= lower) | (time >= upper))
        return(matrix(0, nrow=1, ncol=ncol(base)))

    Exp_gamma <- exp(intensity_coef*prev_outcome)

    pmf <- pmf_estimator(Xb_ind)

    # Everything above this does not depend on alpha and could be optimized
    E_exp_alphaY <- crossprod( exp(alpha*y), pmf )
    E_Yexp_alphaY <- crossprod( y*exp(alpha*y), pmf )
    E_Y_past <- E_Yexp_alphaY/E_exp_alphaY

    Term1_unweighted <-
        (outcome-E_Y_past)/
        (baseline_lambda*Exp_gamma* exp(-alpha*outcome)*E_exp_alphaY)

    pcoriaccel_evaluate_basis(base, time) * c(Term1_unweighted)
}

compute_influence_term_1_for_all <-
    function(
        times_all,
        X_all,
        outcome_all,
        prev_outcome_all,
        baseline_intensity_all,
        alpha,
        intensity_coef,
        outcome_coef,
        base,
        bandwidth,
        kernel
    ){
        if(getOption('PCORI::do_arg_checks', TRUE))
            assert_that(
                is.vector(times_all), is.double(times_all), all(times_all > 0),
                is.matrix(X_all), is.double(X_all), nrow(X_all) == length(times_all),
                is.vector(outcome_all), is.double(outcome_all), length(outcome_all) == nrow(X_all),
                is.vector(prev_outcome_all), is.double(prev_outcome_all), length(prev_outcome_all) == nrow(X_all),
                is.vector(baseline_intensity_all), is.double(baseline_intensity_all), length(baseline_intensity_all) == nrow(X_all),
                rlang::is_scalar_double(intensity_coef),
                is.vector(outcome_coef), is.double(outcome_coef), length(outcome_coef) == ncol(X_all),
                is(base, "SplineBasis"),
                rlang::is_scalar_double(bandwidth), bandwidth > 0,
                rlang::is_string(kernel)
            )


        y_seq <- sort(unique(outcome_all))
        Xbeta <- X_all %*% outcome_coef
        pmf_estimator <-
            function(x){
                pcoriaccel_estimate_pmf(
                    Xb = Xbeta, Y = outcome_all,
                    xi = x, y_seq = y_seq,
                    h = bandwidth,
                    kernel = kernel
                )
            }


        term1 <- matrix(NA_real_, nrow=nrow(X_all), ncol=dim(base)[2])
        for (i in seq_len(nrow(X_all))) {
            term1[i,] <-
            compute_influence_term_1_at_timepoint(
                time = times_all[i],
                outcome = outcome_all[i],
                prev_outcome = prev_outcome_all[i],
                Xb_ind = Xbeta[i],
                pmf_estimator = pmf_estimator,
                intensity_coef = intensity_coef,
                alpha = alpha,
                baseline_lambda = baseline_intensity_all[i],
                base = base,
                y = y_seq
            )
        }
        return(term1)
    }


compute_influence_term_2_for_individual <-
    function(
        times_ind,
        X_ind,
        X_all,
        Y_all,
        slope,
        alpha,
        beta,
        base,
        bandwidth,
        tol = 1e-6
        ){

        if (getOption('PCORI::do_arg_checks', TRUE))
            assert_that(
                is.numeric(times_ind),
                is.matrix(X_ind), is.double(X_ind), length(times_ind) == nrow(X_ind),
                is.matrix(X_all), is.double(X_all), ncol(X_all) == ncol(X_ind),
                is.vector(Y_all), is.double(Y_all), length(Y_all) == nrow(X_all),
                is.vector(slope), is.double(slope), ncol(X_ind) == length(slope),
                rlang::is_scalar_double(alpha),
                is.vector(beta), is.double(beta), length(beta) == ncol(X_all),
                is(base, "SplineBasis"),
                rlang::is_scalar_double(bandwidth), bandwidth > 0,
                rlang::is_scalar_double(tol), tol > 0,
                all(times_ind >= 0)
            )

        pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix(
            X = X_all,
            Y = Y_all,
            individual_X = X_ind,
            times = times_ind,
            x_slope = slope,
            alpha = alpha,
            beta = beta,
            spline_basis = base,
            bandwidth = bandwidth,
            tol = tol)

}
compute_influence_term_2_for_all_patients <-
    function(
        outcome.model,
        X_all, Y_all, outcome_coef,
        data_all,
        alpha,
        base,
        tol = 1e-6
    ){
        integration_data <- rlang::inject(
            model.frame(
                terms(outcome.model),
                data = data_all |>
                    dplyr::select(-..prev_outcome..) |>
                    dplyr::mutate(
                        ..prev_outcome..  = ..outcome..,
                        ..delta_time..    = 0
                    )
                ,
                id   = ..id..,
                time = ..time..
            )
        )
        integration_X <- model.matrix( outcome.model, data=integration_data)
        time <- integration_data[["(time)"]]
        ids  <- integration_data[["(id)"  ]]
        uids <- unique(ids)


        slope <- model.matrix(
            outcome.model,
            data=tibble(
                ..outcome..=0,
                ..time..=c(0, 1),
                ..prev_outcome..=0,
                ..delta_time..=c(0,1)
            )
        ) |> apply(2, diff)

        term2 <- matrix(NA_real_, nrow=length(uids), ncol=dim(base)[2])
        fcnt <- vector('list', length(uids))
        estim.prec <- vector('list', length(uids))
        for (i in seq_along(uids)) {
            result <- compute_influence_term_2_for_individual(
                times_ind = time[ids == uids[i]],
                X_ind = integration_X[ids == uids[i],,drop=FALSE],
                X_all = X_all,
                Y_all = Y_all,
                slope = slope,
                alpha = alpha,
                beta = outcome_coef,
                base = base,
                bandwidth = outcome.model$bandwidth,
                tol = tol
            )

            term2[i,] <- result
            fcnt[[i]] <- attr(result, 'fcnt')
            estim.prec[[i]] <- attr(result, 'estim.prec')
        }
        structure(term2, id = uids, fcnt = fcnt, estim.prec = estim.prec)
    }

compute_influence_terms <-
    function(
        data_all_with_transforms, #< vector of times for all observations
        base, # Spline basis
        alpha, # Sensitivity, singular alpha value
        outcome.model, # outcome model
        intensity_coef, # Coefficient(s) from the intensity model
        baseline_intensity_all, # Baseline intensity for all patients
        tol = 1e-6
    ){
        followup.data <- filter(data_all_with_transforms, ..time.. > 0)
        X_all <- model.matrix(outcome.model)
        Y_all <- model.response(model.frame(outcome.model))
        id <- pull(filter(data_all_with_transforms, ..time.. > 0), ..id..)
        ids <- sort(unique(pull(data_all_with_transforms, ..id..)))
        outcome_coef <- outcome.model$coef


        term1 <-
            compute_influence_term_1_for_all(
                X_all = model.matrix(outcome.model, data = followup.data),
                times_all = pull(followup.data, ..time..),
                outcome_all = pull(followup.data, ..outcome..),
                prev_outcome_all = pull(followup.data, ..prev_outcome..),
                baseline_intensity_all = baseline_intensity_all,
                alpha = alpha,
                intensity_coef = intensity_coef,
                outcome_coef = outcome_coef,
                base = base,
                bandwidth = outcome.model$bandwidth,
                kernel = attr(outcome.model, 'kernel')
            )

        term2 <- compute_influence_term_2_for_all_patients(
            outcome.model = outcome.model,
            X_all=X_all, Y_all=Y_all, outcome_coef=outcome_coef,
            data_all = data_all_with_transforms,
            alpha = alpha,
            base = base,
            tol=tol
        )
        list(
            id = ids,
            alpha = alpha,
            term1 = map(ids, \(.){
                colSums(term1[id == .,, drop=FALSE])
            }) |> do.call(rbind, args=_),
            term2=term2
        )
    }

