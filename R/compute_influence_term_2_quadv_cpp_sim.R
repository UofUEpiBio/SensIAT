compute_influence_term_2_quadv_cpp_sim <-
    function(
        df_i,
        alpha,
        outcome.model,
        base,
        variables,
        centering,
        ...,
        tol=.Machine$double.eps^(1/4)
    ){
        assert_that(
            is.data.frame(df_i),
            is(base, 'SplineBasis')
        )

        # a <- base@knots[[base@order]]
        # b <- tail(base@knots, base@order)[[1]]
        #
        # times <- pull(df_i, variables$time)
        # times <- unique(c(a, times[a < times & times < b], b))
        #
        # mf <- model.frame(outcome.model)
        # Xi <- model.matrix(terms(outcome.model), data = mf)
        # beta <- outcome.model$coef
        # Yi <- model.response(mf)
        # y <- sort(unique(Yi))
        # Xb <- Xi %*% beta
        #
        # period.integrals <-
        #     purrr::map2(head(times, -1), tail(times, -1), \(lower,upper){
        #         lower_df <- impute_patient_df(lower, df_i,
        #                                       variables = variables,
        #                                       centering = centering,
        #                                       right = FALSE)
        #         upper_df <- impute_patient_df(upper, df_i,
        #                                       variables=variables,
        #                                       centering = centering,
        #                                       right = TRUE)
        #
        #         xb_lower <- model.matrix(terms(outcome.model), data=lower_df) %*% outcome.model$coef
        #         xb_upper <- model.matrix(terms(outcome.model), data=upper_df) %*% outcome.model$coef
        #
        #         pcoriaccel_integrate_simp(\(time){
        #             a <- (time - lower)/(upper-lower)
        #
        #             xb_time <- (1-a)*xb_lower + a*xb_upper
        #
        #             pmf <- pcoriaccel_estimate_pmf(Xb, Yi, xb_time, y, outcome.model$bandwidth)
        #
        #             E_exp_alphaY  <- sum(   exp(alpha*y)*pmf )
        #             E_Yexp_alphaY <- sum( y*exp(alpha*y)*pmf )
        #             ev <- E_Yexp_alphaY/E_exp_alphaY
        #
        #             B <- evaluate(base, time)
        #             return(B*ev)
        #         }, lower, upper, tol=tol)
        #     })
        #
        #
        # map(period.integrals, getElement, 'Q') |>
        #     reduce(`+`) |>
        #     structure(
        #         fcnt = purrr::map_int(period.integrals, getElement, 'fcnt'),
        #         estim.prec = purrr::map_dbl(period.integrals, getElement, 'estim.prec')
        #     )




    }
