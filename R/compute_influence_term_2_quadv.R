compute_influence_term_2_quadv <-
function(
    df_i,
    expected_value,
    base,
    variables,
    ...,
    tol=.Machine$double.eps^(1/4)
){
    assert_that(
        is.data.frame(df_i),
        is.function(expected_value),
        is(base, 'SplineBasis')
    )

    a <- base@knots[[base@order]]
    b <- tail(base@knots, base@order)[[1]]

    patient.df <- df_i |>
        filter(
            !!a <= !!variables$time,
            !!variables$time <= !!b
        )

    times <- unique(c(a, pull(patient.df, variables$time), b))
    periods <-
        tibble(
            lower=head(times, -1),
            upper=tail(times, -1)
        )

    period.integrals <-
        pmap(periods, \(lower,upper){
            pracma::quadv(\(time){
                imputed.data <-
                    if_else(
                        time == lower,
                        impute_patient_df(time, df_i, variables=variables, ..., right = FALSE),
                        impute_patient_df(time, df_i, variables=variables, ..., right = TRUE)
                    )
                ev <- expected_value(imputed.data)[,]
                B <- evaluate(base, time)
                return(B*ev)
            }, lower, upper, tol=tol)
        })


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
     fit_SensIAT_within_group_model(
         group.data = filter(ARC_data, Trt=='home_visits'),
         outcome_modeler = SensIAT_sim_outcome_modeler,
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
