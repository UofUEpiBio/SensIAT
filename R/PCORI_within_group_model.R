#' Produce fitted model for group (treatment or control)
#'
#' Produces a fitted model that may be used to produce estimates of mean and
#' variance for the given group.
#'
#' This function should be agnostic to whether it is being provided a
#' treatment or control group.
#'
#' @param group.data The data for the group that is being analyzed.
#'          Preferrably passed in as a single tibble that internally is
#'          subsetted/filtered as needed.
#' @param End The end time for this data analysis, we need to set the default value as the
#'           max value of the time
#' @param kernel
#' @param nmk
#' @param outcome_modeler A separate function that may be swapped out to switch
#'          between negative-binomial, single index model, or another we will
#'          dream up in the future.
#' @param ... add parameters as needed or use this to pass forward into the
#'          outcome_modeler.
#'
#' @return
#'      Should return everything needed to define the fit of the model.
#'      This can then be used for producing the estimates of mean, variance,
#'      and in turn treatment effect.
#'
#' @export
#'
#' @examples
#'
#' fitted.trt.nb <-
#'     fit_PCORI_within_group_model(
#'         group.data = filter(ARC_data, Trt=='home_visits'),
#'         id.var = elig_pid,
#'         outcome.var = Asthma_control,
#'         time.var = time,
#'         intensity.bandwidth = 30,
#'         End = 830,
#'         outcome_modeler = glm.nb
#'     )
#'
#' fitted.trt <-
#'     fit_PCORI_within_group_model(
#'         group.data = filter(ARC_data, Trt=='home_visits'),
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         id.var = elig_pid,
#'         outcome.var = Asthma_control,
#'         time.var = time,
#'         intensity.bandwidth = 30,
#'         End = 830
#'     )
#'
fit_PCORI_within_group_model <- function(
        group.data,
        outcome_modeler,
        id.var,
        outcome.var,
        time.var,
        intensity.covariates = ~1,
        intensity.bandwidth,
        outcome.covariates = ~-1,
        End = max({{time.var}}, na.rm = TRUE) + 1,
        ...
){
    id.var <- rlang::ensym(id.var)
    outcome.var <- rlang::ensym(outcome.var)
    time.var <- rlang::ensym(time.var)
    outcome_modeler <- match.fun(outcome_modeler)
    End <- rlang::eval_tidy({{End}}, data = group.data, env =parent.frame())

    # For example, if we want to fit models for the treatment group
    # Then the group.data should includes parameter "End"
    # Dataframe: "data_baseline_hv", "data_visits_hv", "data_survival_hv"

    group.data2 <- filter(group.data, time <= End)

    data_formatted <- formatting_fn(df = group.data2,
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
    # data_baseline_hv <- data_formatted[[1]]
    data_visits_hv <- data_formatted[[2]]
    # data_survival_hv <- data_formatted[[3]]

    u_hv <- group.data2 |> select(!!id.var) |> dplyr::distinct() |> pull()
    N <- pull(dplyr::summarize(group.data2, n_distinct(!!id.var)))
    # N <- dim(data_baseline_hv)[1]

    ###########################   Model 1: Intensity model  ######################

    # Andrew [@halpo]: the user can only change this part, "Prev_outcome+strata(Visit_number)"

    ######   Andersen-Gill model stratifying by assessment number

    model.data <-
        rlang::inject(!!outcome.var ~ !!id.var + !!time.var + !!rlang::f_rhs(intensity.covariates)) |>
        model.frame(data=dplyr::filter(group.data, (!!time.var) < !!End)) |>
        dplyr::arrange(!!id.var, !!time.var) |>
        dplyr::group_by(!!id.var) |>
        dplyr::mutate(
            visit.number = seq_along(!!time.var)
        ) |>
        dplyr::ungroup() |>
        tidyr::complete(!!id.var, visit.number, fill = tibble::lst(!!time.var := !!End)) |>
        dplyr::group_by(!!id.var) |>
        dplyr::arrange(!!id.var, visit.number) |>
        dplyr::mutate(
            "lag({outcome.var})" := dplyr::lag(!!outcome.var, order_by = !!time.var),
            "lag({time.var})"    := dplyr::lag(!!time.var, order_by =  !!time.var),
            dplyr::across(time.var, dplyr::coalesce, !!End)
        ) |>
        filter(!(visit.number > 1 & is.na(dplyr::lag(!!outcome.var))))

    intensity.model <-
        rlang::inject(coxph(
            Surv(!!rlang::sym(glue::glue("lag({time.var})")),
                 !!time.var,
                 !is.na(!!outcome.var))~(!!rlang::sym(glue::glue("lag({outcome.var})")))+strata(visit.number),
            id = !!id.var,
            data = filter(model.data, !!time.var > 0)
        ))

    # gamma <- Int_model$coefficients # parameter lambda in lambda(t, O(t))
    # we need this gamma value


    ############  Estimated baseline intensities for each stratum ############
    {


        # data_surv <- survfit(intensity.model, newdata=data.frame(Prev_outcome=0))
        data_surv <- survfit(intensity.model, newdata=tibble("lag({outcome.var})" := 0))
        strata <- data_surv$strata

        # Andrew [@halpo]: change this part for a general version (use that general v)

        # the following is to find the baseline intensity function for each strata k

        cumhaz.data <-
            tibble(time = data_surv$time,
                   cumhaz = data_surv$cumhaz,
                   strata = factor(rep(names(strata), strata), levels = names(strata))
                   ) |>
            dplyr::group_by(strata) |>
            dplyr::mutate(hazard = cumhaz - dplyr::lag(cumhaz, default=0, order_by = time))
        base_intens <- cumhaz.data |>
            dplyr::group_by(strata) |>
            dplyr::group_map(rlang::inject(~purrr::map_dbl(seq.int(!!End), lambda0_fn, b=intensity.bandwidth, surv = .x)))


        tmp <- outer(data_surv$time, seq.int(End), `-`)

        base_intens <- apply(
            0.75*(1 - (tmp/intensity.bandwidth)**2) * (abs(tmp) < intensity.bandwidth) * cumhaz.data$hazard,
            2, tapply, cumhaz.data$strata, sum
        )

        # v1=strata[1]
        # cumhaz_v1=data.frame(time <- data_surv$time[1:v1],
        #                      cumhaz <- data_surv$cumhaz[1:v1])
        # base_intens_v1=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)
        #
        # v2=strata[2]
        # cumhaz_v2=data.frame(time=data_surv$time[(v1+1):(v1+v2)],
        #                      cumhaz=data_surv$cumhaz[(v1+1):(v1+v2)])
        # base_intens_v2=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)
        #
        # v3=strata[3]
        # cumhaz_v3=data.frame(time=data_surv$time[(v1+v2+1):(v1+v2+v3)],
        #                      cumhaz=data_surv$cumhaz[(v1+v2+1):(v1+v2+v3)])
        # base_intens_v3=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)
        #
        # v4=strata[4]
        # cumhaz_v4=data.frame(time=data_surv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
        #                      cumhaz=data_surv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
        # base_intens_v4=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)
    }


    #############  Model 2: Outcome model - single index model   ############
    {

        # Andrew [@halpo]: keep the single index model as a flexible model

        # Question: Does this make sense to scale the Lag_time?  should we be lagging scaled time instead?
        outcome.formula <- rlang::inject(
            !!outcome.var~
                ns(!!rlang::sym(glue::glue("lag({outcome.var})")), df=3) +
                scale(!!time.var) +
                scale(!!time.var - !!rlang::sym(glue::glue("lag({time.var})"))) +
                !!rlang::f_rhs(outcome.covariates)
        )
        # data_visits_hv1 <- cbind(data_visits_hv,
        #                          ns1 = ns(data_visits_hv$Prev_outcome, df = 3)[, 1],
        #                          ns2 = ns(data_visits_hv$Prev_outcome, df = 3)[, 2],
        #                          ns3 = ns(data_visits_hv$Prev_outcome, df = 3)[, 3],
        #                          time_scale = scale(data_visits_hv$time),
        #                          Lag_time_scale = scale(data_visits_hv$Lag_time))

        # time_mean <- mean(data_visits_hv1$time)
        # time_sd <- sd(data_visits_hv1$time)
        # Lag_time_mean <- mean(data_visits_hv1$Lag_time)
        # Lag_time_sd <- sd(data_visits_hv1$Lag_time)

        centering.statistics <-
        dplyr::summarize( dplyr::ungroup(dplyr::filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))
                        , "mean({time.var})" := mean(!!time.var)
                        , "sd({time.var})" := sd(!!time.var)
                        , "mean({time.var} - lag({time.var}))" := mean(!!time.var - !!rlang::sym(glue::glue("lag({time.var})")))
                        , "sd({time.var} - lag({time.var}))" := sd(!!time.var - !!rlang::sym(glue::glue("lag({time.var})")))
                        )


        # Xi <- data.frame(data_visits_hv1$ns1,
        #                  data_visits_hv1$ns2,
        #                  data_visits_hv1$ns3,
        #                  data_visits_hv1$time_scale,
        #                  data_visits_hv1$Lag_time_scale)
        #
        # Xi <- as.matrix(Xi)
        # Yi <- as.matrix(data_visits_hv1$Asthma_control, ncol = 1)

        Xi <- model.matrix(outcome.formula, data = dplyr::filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))
        Yi <- model.response(model.frame(outcome.formula, data = dplyr::filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))) |>
                as.matrix(ncol=1)

        # Andrew [@halpo]: add these functions into the R package
        # new method to get the SIM's estimation
        # SDR1 <- cumuSIR_new(X = Xi, Y = Yi)
        # Outcome_model <- SIDRnew(X = Xi, Y = Yi, initial = SDR1$basis[, 1], kernel = "dnorm", method = "nmk")
    }
    outcome.model <- outcome_modeler(outcome.formula, data = dplyr::filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))

    structure(list(
        intensity_model = intensity.model,
        baseline_intensity = base_intens,
        outcome_model = outcome.model,
        outcome_model_centering = centering.statistics,
        data = model.data,
        variables = list(
            id = id.var,
            time = time.var,
            outcome = outcome.var
        ),
        End = End
    ), class = "PCORI_within_group_model")
}



#' Prediction method for `PCORI_within_group_model` objects
#'
#' @param object a `PCORI_within_group_model` object.
#' @param time The time(s) to evaluate at.
#' @param alpha The alpha values to evaluate at.
#'
#'  Evaluate the fitted model, `object`, at each combination of `time` and
#'  `alpha`.
#'
#'
#' @return A tibble/data.frame with the components: `time`, `alpha`, `mean`, `var`.
#' where `time` and `alpha` are the combinations of the respective input(s) and
#' `mean` and `var` are the estimated mean and variance of the response for the
#' given model.
#'
#'
#' @export
#'
#' @examples
#'
#' fitted.trt.nb <-
#'     fit_PCORI_within_group_model(
#'         group.data = filter(ARC_data, Trt=='home_visits'),
#'         id.var = elig_pid,
#'         outcome.var = Asthma_control,
#'         time.var = time,
#'         intensity.bandwidth = 30,
#'         End = 830,
#'         outcome_modeler = glm.nb
#'     )
#' predict(fitted.trt.nb, time = c(90, 180), alpha = c(-0.6, -0.3, 0, 0.3, 0.6))
#'
`predict.PCORI_within_group_model` <-
function(object, time, alpha,
         ){
    #TODO

    spline_seq =


    # Input parameter:
    # spline_seq: the whole time interval End
    # p: dimension of parameter beta in E[Y(t)]'s model
    # B_t_matrix, V, V_inverse, Weights_term2 can be calculated by another function

    # Yujing : check this, see whether we can change this for some existing function - need to do
    B_t_matrix=matrix(sapply(spline_seq,spline_fn),byrow=FALSE, nrow=p)

    ####  we approximated the integral V=\int_t B(t)B(t)'dt using sums of rectangles of width 1
    V=B_t_matrix%*%t(B_t_matrix)
    V_inverse=solve(V)

    Weights_term2=V_inverse%*%B_t_matrix


    #######   Get each subject's baseline intensity at each of their own visit times  #######
    ####  assessments in the inference period [a,b]
    Visits_df <- filter(data_visits_hv1, time >= min(spline_seq),time <= max(spline_seq))
    K <- dim(Visits_df)[1]

    baseline_lambda <- rep(NA,K)

    for(k in 1:K){

        visit_number=Visits_df$Visit_number[k]
        time=Visits_df$time[k]


        if(visit_number==1){

            baseline_lambda[k]=base_intens_v1[time]
        }

        if(visit_number==2){

            baseline_lambda[k]=base_intens_v2[time]
        }

        if(visit_number==3){

            baseline_lambda[k]=base_intens_v3[time]
        }

        if(visit_number==4){

            baseline_lambda[k]=base_intens_v4[time]
        }

    }

    Visits_df <- mutate(Visits_df, baseline_lambda)

    #######   For the given alpha  #######

    Visits_df_a <- Visits_df

    E_Y_past <- rep(NA,K)
    E_exp_alphaY <- rep(NA,K)
    Exp_gamma <- rep(NA,K)

    # Yujing: send  Cond_mean_fn_single2 to Andrew - done
    ############## This part is for single index model ##############
    for(k in 1:K){
        df_k <- Visits_df_a[k, ]
        Exp_gamma[k] <- exp(gamma*df_k$Prev_outcome)

        temp <- Cond_mean_fn_single2(alpha,
                                     X = Xi,
                                     Y = Yi,
                                     x = c(df_k$ns1,
                                           df_k$ns2,
                                           df_k$ns3,
                                           df_k$time_scale,
                                           df_k$Lag_time_scale),
                                     beta = outcome_model$coef,
                                     bandwidth = outcome_model$bandwidth)

        E_Y_past[k] <- temp[[1]]
        E_exp_alphaY[k] <- temp[[2]]
    }

    Visits_df_a <- mutate(Visits_df_a, Exp_gamma, E_Y_past, E_exp_alphaY)

    Visits_df_a <- mutate(Visits_df_a, Term1_unweighted=(Asthma_control-E_Y_past)/
                              (baseline_lambda*Exp_gamma*exp(-alpha*Asthma_control)*E_exp_alphaY) )

    ##########  Term 1 of the influence function:
    ##########   the kth column is the term corresponding to visit k, in Term 1
    Term1_mat <- matrix(nrow=p, ncol=K)

    for(k in 1:K){
        time_k <- Visits_df_a$time[k]
        spline_k <- matrix(spline_fn(time_k),ncol=1)
        Term1_mat[,k] <- (V_inverse%*%spline_k) * Visits_df_a$Term1_unweighted[k]
    }


    ############  Term 2 of the influence function:
    ############   the ith column is Term 2 for participant i

    Term2_mat=matrix(nrow=p, ncol=N_hv)

    # this part matches with the outcome model - single index model
    data_survival_hv1 <- cbind(data_survival_hv,
                               ns1 = ns(data_survival_hv$Prev_outcome, df = 3)[, 1],
                               ns2 = ns(data_survival_hv$Prev_outcome, df = 3)[, 2],
                               ns3 = ns(data_survival_hv$Prev_outcome, df = 3)[, 3],
                               time_scale = (data_survival_hv$time - time_mean) / time_sd,
                               Lag_time_scale = (data_survival_hv$Lag_time - Lag_time_mean) / Lag_time_sd)

    for(i in 1:N_hv){
        print(i)
        df_i1=filter(data_survival_hv1, elig_pid==u_hv[i])

        ############## change this part for single index model ##############
        {
            # for the time interval spline_seq, generate the covariate matrix X
            length_time <- length(spline_seq)
            df_est <- data.frame(time = spline_seq,
                                 Lag_time = rep(0, length_time),
                                 ns1 = rep(0, length_time),
                                 ns2 = rep(0, length_time),
                                 ns3 = rep(0, length_time),
                                 time_scale = rep(0, length_time),
                                 Lag_time_scale = rep(0, length_time))

            # new version
            {
                for(k in 1:length_time){
                    t <- spline_seq[k]

                    min_time <- min(df_i1$time)
                    max_time <- max(df_i1$time)

                    if(t <= min_time){
                        df_est[k, 2:7] <- c(t,
                                            # df_i1$Prev_outcome[1],
                                            df_i1$ns1[1],
                                            df_i1$ns2[1],
                                            df_i1$ns3[1],
                                            (t - time_mean) / time_sd,
                                            (t - Lag_time_mean) / Lag_time_sd)
                    }else if(t >= max_time){
                        temp <- df_i1$Asthma_control[nrow(df_i1)]
                        df_est[k, 2:7] <- c(t - max_time,
                                            # temp,
                                            ns(temp + 0.01 , df = 3)[1],
                                            ns(temp + 0.01, df = 3)[2],
                                            ns(temp + 0.01, df = 3)[3],
                                            (t - time_mean) / time_sd,
                                            (t - max_time - Lag_time_mean) / Lag_time_sd)
                    }else{
                        # t within min_time and max_time
                        t_index1 <- max(which(df_i1$time < t))
                        t_index2 <- min(which(df_i1$time > t))
                        df_est[k, 2:7] <- c(t - df_i1$time[t_index1],
                                            # df_i1$Prev_outcome[t_index2],
                                            df_i1$ns1[t_index2],
                                            df_i1$ns2[t_index2],
                                            df_i1$ns3[t_index2],
                                            (t - time_mean) / time_sd,
                                            (t - df_i1$time[t_index1] - Lag_time_mean) / Lag_time_sd)
                    }

                }
            }

            # start <- Sys.time()
            Time_means_single  <- matrix(0, nrow = length_time, ncol = 1)

            start <- Sys.time()
            for(l in 1:length_time){
                temp <- Cond_mean_fn_single2(alpha,
                                             X = Xi,
                                             Y = Yi,
                                             x = c(df_est$ns1[l],
                                                   df_est$ns2[l],
                                                   df_est$ns3[l],
                                                   df_est$time_scale[l],
                                                   df_est$Lag_time_scale[l]),
                                             beta = outcome_model$coef,
                                             bandwidth = outcome_model$bandwidth)
                Time_means_single[l] <- temp[[1]]
            }
            # end <- Sys.time()
            # end - start
            Term2_mat[,i] <- Weights_term2%*%Time_means_single
        }
    }

    #######   Subject-specific p x 1 influence function  IF(O_i)
    IF_mat <- matrix(nrow=p, ncol=N_hv)

    for(i in 1:N_hv){

        # w=which(Visits_df_a$elig_pid==u[i]) # the previous code
        w = which(Visits_df_a$elig_pid == u_hv[i]) # the new code

        if(length(w) ==0){ temp = 0 }
        if(length(w) ==1){ temp=Term1_mat[,w] }
        if(length(w) > 1){ temp=rowSums(Term1_mat[,w]) }

        IF_mat[,i]=temp + Term2_mat[,i]

    }

    ########  Target parameter (p x 1)
    Beta_hat <- rowMeans(IF_mat)
    # estimation of beta in E[Y(t)] = s[beta %*% B(t)] # s(\cdot): identity link

    ########  IF-based variance estimation
    Var_array <- array(dim=c(p, p, N_hv))

    for(i in 1:N_hv){
        temp <- IF_mat[ , i] - Beta_hat
        Var_array[ , , i] <-  temp%*%t(temp)
    }
    Var_beta <- (1/N_hv^2)*rowSums(Var_array, dims=2) # estimation of variance for beta


    ####### The previous part is based on the "alpha" parameter
    ####### The following part is based on the "time" parameter

    #######  Target-time means and IF-based variance estimates

    #######   For the given time  #######
    B_t <- matrix(spline_fn(time), ncol = p)

    mean_t  <- B_t%*%Beta_hat # estimation of E[Y(t); alpha]
    Var_t  <- B_t%*%Var_beta%*%t(B_t) # estimation of Var(Y(t); alpha)

    return(data.frame(alpha, time, mean_t, Var_t))
}
