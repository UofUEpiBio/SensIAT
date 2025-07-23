compute_sim_influence_term_1_at_timepoint <- function(
        outcome,
        Xb_ind, #< individual-level linear predictor to be passed to pmf_estimator
        pmf_estimator,
        alpha, #< sensitivity parameter.
        base,
        intensity_weight,
        y = sort(unique(outcome))
){
    if (getOption('SensIAT::do_arg_checks', TRUE))
        assert_that(
            rlang::is_scalar_double(outcome),
            rlang::is_scalar_double(Xb_ind),
            is.function(pmf_estimator),
            rlang::is_scalar_double(alpha)
        )

    pmf <- pmf_estimator(Xb_ind)

    # Everything above this does not depend on alpha and could be optimized
    E_exp_alphaY  <- as.vector(pmf %*%    exp(tcrossprod(y, alpha)) )
    E_Yexp_alphaY <- as.vector(pmf %*% (y*exp(tcrossprod(y, alpha))))

    # E_exp_alphaY <- crossprod( exp(alpha*y), pmf )[,,drop=TRUE]
    # E_Yexp_alphaY <- crossprod( y*exp(alpha*y), pmf )[,,drop=TRUE]
    E_Y_past <- E_Yexp_alphaY/E_exp_alphaY

    (outcome-E_Y_past)/
        (intensity_weight * exp(-alpha*outcome)*E_exp_alphaY)
}


compute_sim_influence_term_1_for_all <-
    function(
        times_all,
        X_all,
        outcome_all,
        alpha,
        outcome_coef,
        intensity_weights,
        base,
        bandwidth,
        kernel
    ){
        if(getOption('SensIAT::do_arg_checks', TRUE))
            assert_that(
                is.numeric(times_all),
                is.matrix(X_all), is.double(X_all),
                length(times_all) == nrow(X_all),
                is.vector(outcome_all), is.double(outcome_all),
                length(outcome_all) == nrow(X_all),
                is.vector(intensity_weights), is.double(intensity_weights),
                length(intensity_weights) == nrow(X_all),
                is.vector(outcome_coef), is.double(outcome_coef),
                length(outcome_coef) == ncol(X_all),
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

        lower <- base@knots[base@order]
        upper <- base@knots[length(base@knots) - base@order+1]

        term1 <- matrix(NA_real_, nrow=nrow(X_all), ncol=dim(base)[2])
        for (i in seq_len(nrow(X_all))) {
            if((times_all[i] <= lower) || (times_all[i] >= upper)){
                term1[i,] <- 0
                next
            }

            term1[i,] <- pcoriaccel_evaluate_basis(base, times_all[i]) *
                compute_sim_influence_term_1_at_timepoint(
                    outcome = outcome_all[i],
                    Xb_ind = Xbeta[i],
                    pmf_estimator = pmf_estimator,
                    alpha = alpha,
                    base = base,
                    y = y_seq,
                    intensity_weight = intensity_weights[i]
                )
        }
        return(term1)
    }


compute_sim_influence_term_2_for_individual <-
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
                compute_sim_influence_term_2_for_one_patient_fixed_approximation(
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

compute_sim_influence_term_2_for_one_patient_fixed_approximation <-
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



compute_sim_influence_term_2_for_all_patients <-
    function(
        outcome.model,
        integration_data, #< includes baseline with prev_outcome, and prev_time
                          #< replaced with the current values.
        alpha,
        base,
        id,time, #< lazy evaluation
        time.vars = character(), #< Not lazy evaluation, additional time variables
        ...
    ){
        id <- rlang::ensym(id)
        time <- rlang::ensym(time)
        terms <- terms(outcome.model)
        mf <- rlang::inject(model.frame(terms, data=integration_data, id = !!id, time=!!time))
        integration_X <- model.matrix( terms, data=mf)
        times <- mf[["(time)"]]
        ids  <- mf[["(id)"]]
        uids <- unique(ids)

        X_all <- model.matrix(outcome.model)
        Y_all <- model.response(model.frame(outcome.model))

        slope <-
            compute_slope(
                outcome.model=outcome.model,
                time.vars = c(rlang::as_string(time), time.vars),
                ... #< IGNORED
            )

        term2 <- matrix(NA_real_, nrow=length(uids), ncol=dim(base)[2])
        fcnt <- vector('list', length(uids))
        estim.prec <- vector('list', length(uids))
        for (i in seq_along(uids)) {
            result <- compute_sim_influence_term_2_for_individual(
                times_ind = times[ids == uids[i]],
                X_ind = integration_X[ids == uids[i],,drop=FALSE],
                X_all = X_all,
                Y_all = Y_all,
                slope = slope,
                alpha = alpha,
                beta = coef(outcome.model),
                base = base,
                bandwidth = outcome.model$bandwidth,
                kernel = attr(outcome.model, 'kernel'),
                ...
            )

            term2[i,] <- result
            fcnt[[i]] <- attr(result, 'fcnt')
            estim.prec[[i]] <- attr(result, 'estim.prec')
        }
        structure(term2, id = uids, fcnt = fcnt, estim.prec = estim.prec)
    }

#' @describeIn compute_influence_terms Optimized method for the single index model.
#' @param tolerance Numeric value indicating the tolerance for integration, default is `.Machine$double.eps^(1/3)`.
#' @param na.action Function to handle missing values, default is `na.fail`.
#' @param time Variable indicating the time variable in the data, by Default will be extracted from the intensity model response.
#' @export
`compute_influence_terms.SensIAT::Single-index-outcome-model` <-
    function(
        outcome.model,
        intensity.model,
        alpha,
        data,
        base,
        tolerance = .Machine$double.eps^(1/3),
        na.action = na.fail,
        id = NULL,
        time = NULL,
        ...
    ){
        if(!missing(id))
            id <- rlang::ensym(id)
        if(missing(id))
            id <- attr(outcome.model, 'id')
        if(!missing(time))
            time <- rlang::ensym(time)
        if(missing(time))
            time <- rlang::sym(rlang::f_lhs(terms(intensity.model))[[3]])

        followup_only <- data |>
            filter(0 < !!time) |>
            select( !!id, !!time,
                  , all_of(all.vars(terms(intensity.model)))
                  , all_of(all.vars(terms(outcome.model)))
            ) |>
            na.omit()
        mf_followup    <- rlang::inject(model.frame(terms(outcome.model),
                                                    data=followup_only,
                                                    time = !!time,
                                                    id = !!id))
        X_followup     <- model.matrix(outcome.model, data=mf_followup)
        Y_followup     <- model.response(mf_followup)
        ids_followup   <- pull(mf_followup, "(id)")
        times_followup <- pull(mf_followup, "(time)")
        uids  <- sort(unique(ids_followup))

        intensity <-
            estimate_baseline_intensity(
                intensity.model = intensity.model,
                data = followup_only
            )
        exp_gamma <- predict(intensity.model, newdata = followup_only, type='risk', reference='zero')
        intensity_weights <- intensity$baseline_intensity * exp_gamma


        term1 <- compute_sim_influence_term_1_for_all(
                X_all                  = X_followup,
                times_all              = times_followup,
                outcome_all            = Y_followup,
                alpha = alpha,
                intensity_weights = intensity_weights,
                outcome_coef = outcome.model$coef,
                base = base,
                bandwidth = outcome.model$bandwidth,
                kernel = attr(outcome.model, 'kernel')
            )
        term1.by.id <- map(uids, \(.){
            colSums(term1[ids_followup == .,, drop=FALSE], na.rm = TRUE)
        }) |> do.call(rbind, args=_)

        outcome <- terms(outcome.model)[[2]]
        integration_data <-
                data |>
                    dplyr::select(-any_of('..prev_outcome..')) |>
                    dplyr::mutate(
                        ..prev_outcome..  = !!outcome,
                        ..delta_time..    = 0,
                        ..prev_time..     = !!time,
                    )

        term2 <- compute_sim_influence_term_2_for_all_patients(
            outcome.model = outcome.model,
            integration_data = integration_data,
            alpha = alpha,
            base = base,
            tol=tolerance,
            id = !!id,
            time = !!time,
            ...
        )
        list(
            id = uids,
            alpha = alpha,
            term1 = term1.by.id,
            term2 = term2
        )
    }
