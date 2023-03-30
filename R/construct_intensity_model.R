#' Construct Intensity model
#'
#' @param survival.data The Survival formatted data limited to the treatment group.
#'                      Must include named components `time`, `Event`, `Prev_time`,
#'                      `Prev_outcome`, and `Visit_number`
#' @param id.var The patient identifier variable in `survival_data`
#'
#' @return A list of lists where at the top most level one list for each stratum,
#'      and for eaach stratum a list with components `cumhaz`, the cummulative hazard,
#'      and `base_intens`, the baseline intensity.
#' @export
#'
#' @examples
#' ARC_formatted <- formatting_fn(
#'     df              = ARC_data,
#'     id_var          = "elig_pid",
#'     treatment_var   = "Trt",
#'     outcome_var     = "Asthma_control",
#'     visitnumber_var = "Visit_number",
#'     time_var        = "time",
#'     V               = 5,
#'     knots           = c(59,59,59,59,260,461,461,461,461),
#'     spline_seq      = 60:460,
#'     last_day = 830
#' )
#' survival.data <- ARC_formatted$Survival |> filter(Trt == "home_visits")
#' intensity <- construct_intensity_model(survival_data)
#' glimpse(intensity)
construct_intensity_model <-
function(
    survival.data,
    id.var = names(survival.data)[[1]]
){
    AG_model <-
        coxph(
            Surv(Prev_time,time,Event)~Prev_outcome+strata(Visit_number),
               id = survival.data[[id.var]],
               data = survival.data)
    AG_surv <- survfit(AG_model,newdata=data.frame(Prev_outcome=0))
    strata <- AG_surv$strata

    End <- max(survival.data$time, na.rm = TRUE)


    v1 <- strata[1]
    cumhaz_v1 <- data.frame(time=AG_surv$time[1:v1],
                            cumhaz=AG_surv$cumhaz[1:v1])
    # head(cumhaz_v1)
    base_intens_v1 <- sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)


    v2 <- strata[2]
    cumhaz_v2 <- data.frame(time=AG_surv$time[(v1+1):(v1+v2)],
                            cumhaz=AG_surv$cumhaz[(v1+1):(v1+v2)])
    # head(cumhaz_v2)
    base_intens_v2 <- sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)


    v3 <- strata[3]
    cumhaz_v3 <- data.frame(time=AG_surv$time[(v1+v2+1):(v1+v2+v3)],
                            cumhaz=AG_surv$cumhaz[(v1+v2+1):(v1+v2+v3)])
    head(cumhaz_v3)
    base_intens_v3 <- sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)


    v4 <- strata[4]
    cumhaz_v4 <- data.frame(time=AG_surv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
                            cumhaz=AG_surv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
    head(cumhaz_v4)
    base_intens_v4 <- sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)

    list(
        list(cumhaz = cumhaz_v1, base_intens = base_intens_v1),
        list(cumhaz = cumhaz_v2, base_intens = base_intens_v2),
        list(cumhaz = cumhaz_v3, base_intens = base_intens_v3),
        list(cumhaz = cumhaz_v4, base_intens = base_intens_v4)
    )
}
