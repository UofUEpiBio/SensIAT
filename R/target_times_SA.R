target_times_SA <-
function(
    data,
    treatment,
    patient_id_var,
    intensity.formula,
    spline_seq = attr(data, "spline_seq"),
    Knots = attr(data, 'knots'),
    End = attr(data, 'End')
){
    trt.expr <- rlang::enquo(treatment)

    p=length(Knots)-4
    B_t_matrix=matrix(sapply(spline_seq,spline_fn, Knots),byrow=FALSE,nrow=p)
    V=B_t_matrix%*%t(B_t_matrix)
    V_inverse=solve(V)

    Weights_term2=V_inverse%*%B_t_matrix

    baseline_trt = filter(data$Baseline, !!trt.expr)
    visits_trt = filter(data$Visits, !!trt.expr)
    survival_trt = filter(data$Survival, !!trt.expr)

    N_trt = nrow(baseline_trt)
    u_trt = unique(survival_trt$elig_pid)

    boot_samp <- sample(N_trt,N_trt,replace=TRUE)

    length(boot_samp)

    df_boot <- data.frame(
        boot_id = seq_along(boot_samp),
        elig_pid = boot_samp
    )

    baseline_boot <- inner_join(df_boot, baseline_trt, multiple='all', by=patient_id_var)
    survival_boot <- inner_join(df_boot, survival_trt, multiple='all', by=patient_id_var)
    visits_boot   <- inner_join(df_boot, visits_trt  , multiple='all', by=patient_id_var)


    ###########################    Intensity model  ############################
    AG_model = coxph(intensity.formula, id=elig_pid, data=df_boot)
    gamma = AG_model$coefficients

    ############  Estimated baseline intensities for each stratum  #############
    AG_surv=survfit(AG_model,newdata=data.frame(Prev_outcome=0))
    strata=AG_surv$strata


    base_intens <- dplyr::bind_cols(
        purrr::map2(AG_surv$strata, cumsum(dplyr::lag(AG_surv$strata, 1, 0)),
                    function(stratum, cumprev){
                        cumhaz <- data.frame( time   = AG_surv$time[(1:stratum)+cumprev]
                                              , cumhaz = AG_surv$cumhaz[(1:stratum)+cumprev]
                        )
                        sapply(1:End,lambda0_fn,b=30,surv=cumhaz) #< not sure if b is constant or should be passed as a parameter
                    }))



    ###########################    Outcome model  ###########################
    visits_trt <- filter(df_boot, Event==1)
    visits_trt2 <- data.frame(elig_pid=u_trt[boot_samp]) |>
        left_join(filter(data$Visits, !!trt.expr), by='elig_pid', multiple='all')

    all.equal(visits_trt, visits_trt2)

    Outcome_model <- glm.nb(outcome.formula,data = visits_trt)
    Outcome_model2 <- glm.nb(outcome.formula,data = visits_trt2)
    theta = Outcome_model$theta

    coefs=Outcome_model$coefficients
    a=coefs[5]+coefs[6]


    ##################   Construct influence function for beta  ###############


}
if(FALSE){
    data <- formatting_fn(
            df              = ARC_data,
            id_var          = "elig_pid",
            treatment_var   = "Trt",
            outcome_var     = "Asthma_control",
            visitnumber_var = "Visit_number",
            time_var        = "time",
            V               = 5,
            knots           = c(59,59,59,59,260,461,461,461,461),
            spline_seq      = 60:460,
            last_day = 830
        )
    trt.expr <- rlang::expr(Trt == "home_visits")
    intensity.formula <- Surv(Prev_time,time,Event)~Prev_outcome+strata(Visit_number)
    outcome.formula   <- (6*Asthma_control)~ns(Prev_outcome,df=3) + time + Lag_time
    patient_id_var <- 'elig_pid'


    target_times_SA(
        data,
        Trt == "home_visits",
        "elig_pid",
        intensity.formula
    )
}







