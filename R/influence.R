compute_influence_term_1_unweighted_at_timepoint <- function(
        time,
        outcome,
        Exp_gamma,
        Xb_ind, #< individual-level covariates to be passed to pmf_estimator
        y,pmf,
        alpha, #< sensitivity parameter.
        baseline_lambda
){
    if (getOption('SensIAT::do_arg_checks', TRUE))
        assert_that(
            rlang::is_scalar_double(time), time > 0,
            rlang::is_scalar_double(outcome),
            rlang::is_scalar_double(Exp_gamma),
            rlang::is_scalar_double(Xb_ind),
            is.numeric(pmf), is.numeric(y), length(pmf) == length(y),
            rlang::is_scalar_double(alpha),
            rlang::is_scalar_double(baseline_lambda)
        )

    # Everything above this does not depend on alpha and could be optimized
    E_exp_alphaY <- crossprod( exp(alpha*y), pmf )
    E_Yexp_alphaY <- crossprod( y*exp(alpha*y), pmf )
    E_Y_past <- E_Yexp_alphaY/E_exp_alphaY

    (outcome-E_Y_past)/
    (baseline_lambda*Exp_gamma* exp(-alpha*outcome)*E_exp_alphaY)
}


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
    if (getOption('SensIAT::do_arg_checks', TRUE))
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
    Term1_unweighted <-
    compute_influence_term_1_unweighted_at_timepoint(
        time = time,
        outcome = outcome,
        Exp_gamma = Exp_gamma,
        Xb_ind = Xb_ind,
        y = y,pmf = pmf,
        alpha = alpha,
        baseline_lambda = baseline_lambda
    )
    return(pcoriaccel_evaluate_basis(base, time) * c(Term1_unweighted))
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
        if(getOption('SensIAT::do_arg_checks', TRUE))
            assert_that(
                is.vector(times_all), is.numeric(times_all), all(times_all > 0),
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
        method = c('adaptive', 'fixed'),
        kernel = c("K2_Biweight", "dnorm"),
        tol = 1e-6,
        delta = NULL,
        resolution = NULL,
        fix_discontinuity = TRUE
        ){
        method <- match.arg(method)

        if (getOption('SensIAT::do_arg_checks', TRUE))
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

        if(method == "adaptive"){
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
                tol = tol,
                kernel = kernel)
        } else
        if(method == "fixed"){
            if(!xor(is.null(delta), is.null(resolution)))
                rlang::abort("When method='fixed', either delta or resolution must be provided but not both.")
            if(is.null(delta)){
                assert_that(assertthat::is.count(resolution))
                delta <- diff(range(base@knots))/resolution
            }
            compute_influence_term_2_for_one_patient_fixed_approximation(
                X = X_all,
                Y = Y_all,
                individual_X = X_ind,
                times = times_ind,
                x_slope = slope,
                alpha = alpha,
                beta = beta,
                base = base,
                bandwidth = bandwidth,
                delta = delta,
                kernel = kernel,
                fix_discontinuity=fix_discontinuity
            )


        } else {
            stop("method must be one of 'adaptive' or 'fixed'")
        }

    }

compute_influence_term_2_for_one_patient_fixed_approximation <-
    function(
        X, Y, individual_X, times,
        x_slope,
        alpha,
        beta,
        base,
        bandwidth = bandwidth,
        delta,
        kernel = c("K2_Biweight", "dnorm"),
        fix_discontinuity = TRUE,
        ...
    ){
        kernel <- match.arg(kernel)

        eval.times <- seq(min(base@knots), max(base@knots), by = delta)
        B <- evaluate(base, eval.times)

        uY <- sort(unique(Y))
        Xb <- X %*% beta
        expected_value <- function(time, lhs=FALSE){
            if(lhs){
                period <- max(which(times - time <= 0))
            } else {
                period <- max(which(times - time <  0))
            }
            xi <- individual_X[period, ,drop=FALSE] + (time - times[period]) * x_slope


            pmf <- pcoriaccel_estimate_pmf(
                Xb=Xb, Y = Y,
                xi = xi %*% beta,
                y_seq = uY,
                h= bandwidth, kernel = kernel)
            if(all(pmf ==0)) return(0)
            E_exp_alphaY <- sum( exp(alpha*uY)*pmf )

            E_Yexp_alphaY <- sum( uY*exp(alpha*uY)*pmf )

            return(E_Yexp_alphaY/E_exp_alphaY)
        }


        ipart <- matrix(NA, nrow = length(eval.times)-1, ncol=ncol(B))
        if(!fix_discontinuity)
            ev_right <- expected_value(eval.times[1], lhs=FALSE)
        for(i in seq.int(length(eval.times)-1)){
            if(fix_discontinuity && (eval.times[i] %in% times || i==1)){
                ev_left <- expected_value(eval.times[i], lhs=TRUE)
            } else {
                ev_left <- ev_right
            }
            ev_right <- expected_value(eval.times[i+1], lhs=FALSE)
            ipart[i,] <- (ev_left*B[i,] + ev_right*B[i+1L, ])/2 * delta
        }
        colSums(ipart)
    }


compute_influence_term_2_for_all_patients <-
    function(
        outcome.model,
        X_all, Y_all, outcome_coef,
        data_all,
        alpha,
        base,
        ...
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
                ...
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
        tol = 1e-6,
        ...
    ){
        ..id.. <- ..time.. <- ..prev_outcome.. <- baseline_intensity <- NULL

        followup.data <- filter(data_all_with_transforms, ..time.. > 0) |>
            model.frame(terms(outcome.model), data=_,
                        id = ..id.., time = ..time..,
                        prev_outcome = ..prev_outcome..,
                        intensity = baseline_intensity)

        X_all <- model.matrix(terms(outcome.model), data = followup.data)
        Y_all <- model.response(followup.data)
        id    <- pull(followup.data, '(id)')
        ids   <- sort(unique(pull(data_all_with_transforms, ..id..)))
        outcome_coef <- outcome.model$coef


        term1 <-
            compute_influence_term_1_for_all(
                X_all                  = X_all,
                times_all              = pull(followup.data, "(time)"),
                outcome_all            = Y_all,
                prev_outcome_all       = pull(followup.data, "(prev_outcome)"),
                baseline_intensity_all = pull(followup.data, "(intensity)"),
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
            tol=tol,
            kernel = attr(outcome.model, 'kernel'),
            ...
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

