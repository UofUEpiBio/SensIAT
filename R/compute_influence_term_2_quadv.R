compute_influence_term_2_quadv <- function(df_i, alpha, object, base, ...){
    assert_that(
        is.data.frame(df_i),
        is.numeric(alpha),
        is(object, 'PCORI_within_group_model'),
        is(base, 'SplineBasis')
    )
    expected_value <- \(data, ...){
        matrix(
            object$outcome_model |>
                pcori_conditional_means(..., new.data = data) |>
                pull('E_Y_past'),
            nrow = nrow(data)
        )}
    a <- min(base@knots)
    b <- max(base@knots)
    obase <- orthogonalize(base)

    patient.df <- df_i |>
        filter(
            !!a <= !!object$variables$time,
            !!object$variables$time <= !!b
        )

    periods <-
        tibble(
            lower=c(a, pull(patient.df, object$variables$time)),
            upper=c(pull(patient.df, object$variables$time), b)
        )

    period.integrals <-
        pmap(periods, \(lower,upper, ...){
            pracma::quadv(\(time){
                imputed.data <-
                    if_else(
                        time == lower,
                        impute_patient_df(time, df_i, object, right = FALSE),
                        impute_patient_df(time, df_i, object, right = TRUE)
                    )
                ev <- expected_value(imputed.data, alpha = alpha)[,]
                B <- evaluate(obase, time)
                return(B*ev)
            }, lower, upper, tol=.Machine$double.eps^(1/4))
        })


    map(period.integrals, getElement, 'Q') %>%
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
