compute_influence_term_2_quadv_sim <-
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

    a <- base@knots[[base@order]]
    b <- tail(base@knots, base@order)[[1]]

    times <- pull(df_i, variables$time)
    times <- unique(c(a, times[a < times & times < b], b))

    mf <- model.frame(outcome.model)
    Xi <- model.matrix(terms(outcome.model), data = mf)
    beta <- outcome.model$coef
    Yi <- model.response(mf)
    y <- sort(unique(Yi))
    Xb <- Xi %*% beta

    period.integrals <-
        purrr::map2(head(times, -1), tail(times, -1), \(lower,upper){
            lower_df <- impute_patient_df(lower, df_i,
                                          variables = variables,
                                          centering = centering,
                                          right = FALSE)
            upper_df <- impute_patient_df(upper, df_i,
                                          variables=variables,
                                          centering = centering,
                                          right = TRUE)

            xb_lower <- model.matrix(terms(outcome.model), data=lower_df) %*% outcome.model$coef
            xb_upper <- model.matrix(terms(outcome.model), data=upper_df) %*% outcome.model$coef

            pracma::quadv(\(time, Xb, Yi, xb_lower, xb_upper, y){
                a <- (time - lower)/(upper-lower)

                xb_time <- (1-a)*xb_lower + a*xb_upper

                pmf <- pcoriaccel_estimate_pmf(Xb, Yi, xb_time, y, outcome.model$bandwidth)

                # Fhat <- pcoriaccel_NW(
                #     Xb = Xb, Y = Yi,
                #     xb = xb_time, y = y,
                #     h = outcome.model$bandwidth,
                #     ...)
                # pmf <- diff(c(0, Fhat))
                E_exp_alphaY  <- sum(   exp(alpha*y)*pmf )
                E_Yexp_alphaY <- sum( y*exp(alpha*y)*pmf )
                ev <- E_Yexp_alphaY/E_exp_alphaY

                B <- evaluate(base, time)
                return(B*ev)
            },
            lower, upper,
            Xb, Yi, xb_lower, xb_upper, y,
            tol=tol)
        }
        )


    map(period.integrals, getElement, 'Q') |>
        reduce(`+`) |>
        structure(
            fcnt = purrr::map_int(period.integrals, getElement, 'fcnt'),
            estim.prec = purrr::map_dbl(period.integrals, getElement, 'estim.prec')
        )
}
if(F){
    library(ARCdata)
    fitted.trt.sim <-
     fit_PCORI_within_group_model(
         group.data = filter(ARC_data, Trt=='home_visits'),
         outcome_modeler = PCORI_sim_outcome_modeler,
         id.var = elig_pid,
         outcome.var = Asthma_control,
         time.var = time,
         End = 830
     )
    time.quadv <- system.time({
            pred.quadv <- predict(fitted.trt.sim, time = c(90, 180),
             alpha = c(0),
             intensity.bandwidth = 30,
             knots=c(60,60,60,60,260,460,460,460,460),
             integration.method = 'quadv'
            )
    })
}
