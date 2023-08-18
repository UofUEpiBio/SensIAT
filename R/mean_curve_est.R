




#' Title
#'
#' @param data
#' @param patient_id_var
#' @param treatment
#' @param spline_seq
#' @param Knots
#' @param spline_fn
#' @param ...
#' @param End
#' @param Alpha_seq
#'
#' @return
#' @export
#'
#' @examples
#' estimates <-
#'     estimate_mean_curves(
#'         ARC_data, elig_pid, Trt='home_visits',
#'
#'     )
estimate_mean_curves <-
    function(
    data,
    patient_id_var,
    treatment,
    spline_seq,
    Knots = c( rep(min(spline_seq)-1, 4)
             , mean(range(spline_seq))
             , rep(max(spline_seq)+1, 4)
             ),
    spline_fn = spline_fn,
    ...,
    End = attr(data, 'End'),
    Alpha_seq=c(-0.6,-0.3,0,0.3,0.6)
){
    assertthat::assert_that(
        all(c('Baseline', 'Visits', 'Survival') %in% names(data)),
        msg = "Data is expected to have three components, Baseline, Visits, and Survival."
    )
    patient_id_var <- rlang::ensym(patient_id_var)

    trt.expr <- rlang::enquo(treatment)

    p = length(Knots)-4L #< Length of target parameter beta (the spline parameter)

    #####  the p x length(spline_seq) matrix with columns B(t) defined by `spline_fn`
    B_t_matrix=matrix(sapply(spline_seq, spline_fn, Knots),byrow=FALSE,nrow=p)


    mean_curve_ests(
        baseline = filter(data$Baseline, !!trt.expr),
        survival = filter(data$Survival, !!trt.expr),
        visits   = filter(data$Visits  , !!trt.expr),
        B_t_matrix = B_t_matrix, spline_seq = spline_seq,
        patient_id_var = deparse(patient_id_var),
        End = End,
        Alpha_seq = Alpha_seq,
        knots = Knots,
        ...
    )

}

mean_curve_ests <-
function(
    baseline,
    survival,
    visits,
    patient_id_var,
    B_t_matrix, spline_seq,
    intensity.formula,
    outcome.formula,
    End,
    Alpha_seq,
    Cond_mean_fn,
    knots,
    ...
    ){
    # Nested, intentionally
    estimate_curve_for_one_alpha <-
        function(
        alpha,
        # baseline, survival, visits,
        # patient_id_var,
        # B_t_matrix, spline_seq,
        target = list(
            time = spline_seq,
            B_t = B_t_matrix
        ),
        ...,
        do.var = TRUE
        ) {

        # Visits_df_a=visits

        # E_Y_past=rep(NA,K)
        # E_exp_alphaY=rep(NA,K)
        #
        # Exp_gamma=rep(NA,K)


        # for(k in 1:K){
        #     df_k            = Visits_df_a[k,]
        #     Exp_gamma[k]    = exp(gamma*df_k$Prev_outcome)
        #     mu_k            = predict.glm(Outcome_model,newdata=df_k,type="response")
        #     temp            = Cond_mean_fn(mu=mu_k, theta=Outcome_model$theta, alpha=alpha)
        #     E_Y_past[k]     = temp[[1]]
        #     E_exp_alphaY[k] = temp[[2]]
        # }
        # Visits_df_a=mutate(Visits_df_a,Exp_gamma,E_Y_past,E_exp_alphaY)
        # Visits_df_a=mutate(Visits_df_a,Term1_unweighted=(Asthma_control-E_Y_past)/
        #                        (baseline_lambda*Exp_gamma*exp(-alpha*Asthma_control)*E_exp_alphaY) )

        Exp_gamma = exp(gamma * Visits_df$Prev_outcome)
        mu_k = predict.glm(Outcome_model, newdata = Visits_df, type='response')
        temp = purrr::map(mu_k, Cond_mean_fn, theta = Outcome_model$theta, alpha = alpha)

        names(temp[[1]]) <- c('E_Y_past', 'E_exp_alphaY')

        t2 <- purrr::map_dfr(temp, ~tibble::tibble(E_Y_past = .[[1]], E_exp_alphaY = .[[2]]))

        Visits_df_a <- Visits_df |>
            dplyr::bind_cols(tibble(Exp_gamma), t2) |>
            dplyr::mutate(baseline_lambda = purrr::map2_dbl(time, Visit_number, ~base_intens[[.x, .y]])) |>
            dplyr::mutate(Term1_unweighted=(Asthma_control-E_Y_past)/
                              (baseline_lambda*Exp_gamma*exp(-alpha*Asthma_control)*E_exp_alphaY))



        ##########  Term 1 of the influence function:
        ##########   the kth column is the term corresponding to visit k, in Term 1
        Term1_mat=matrix(nrow=p,ncol=K)

        for(k in 1:K){

            time_k=Visits_df_a$time[k]
            spline_k=matrix(spline_fn(time_k, knots = knots),ncol=1)

            Term1_mat[,k]=(V_inverse%*%spline_k)*Visits_df_a$Term1_unweighted[k]

        }


        ############  Term 2 of the influence function:
        ############   the ith column is Term 2 for participant i

        ##########  Note:  here we leveraged the fact that our outcome model is linear
        ##########    in Time and Lag Time to reduce the number of computations needed
        ##########    in order to get each subject's predicted mean on each day of the
        ##########    period [a,b], so this would need to be modified for different
        ##########    outcome models
        # baseline_trt <- dplyr::filter(data$Baseline, !!trt.expr)
        N_trt = nrow(baseline)
        u_trt = unique(pull(baseline, patient_id_var))

        Term2_mat=matrix(nrow=p,ncol=N_trt)

        for(i in 1:N_trt){

            df_i1=filter(survival,elig_pid==u_trt[i])


            ####  predictors for time = 1
            newdat_i1=data.frame(time=1,Lag_time=1,Prev_outcome=df_i1$Prev_outcome[1])


            df_i2=filter(Visits_df_a,elig_pid==u_trt[i])
            times_i=df_i2$time+1

            if(length(times_i)>0){

                #####  predictors for times just after post-baseline visits in [a,b]
                newdat_i2=data.frame(time=times_i,Lag_time=1,Prev_outcome=df_i2$Asthma_control)

                newdat_i=rbind(newdat_i1,newdat_i2)

            }

            #####  predicted means mu(t*) just after each of participant i's assessments
            mu_i=predict.glm(Outcome_model,newdata=newdat_i,type="response")

            times_i=c(1,times_i)

            time_means_fn=function(t){

                #####  Under the outcome model that we used,
                #####  for each t, mu(t) = mu(t*)exp(a Lag), where t* is the latest
                #####    day-after-an-assessment day on or before day t
                m=max(which(times_i<=t))
                s=t-times_i[m]

                temp=Cond_mean_fn(mu=mu_i[m]*exp(a*s),theta=Outcome_model$theta,alpha=alpha)

                return(temp[[1]])

            }

            Time_means=matrix(sapply(spline_seq,time_means_fn),ncol=1)

            Term2_mat[,i]=Weights_term2%*%Time_means


        }



        #######   Subject-specific p x 1 influence function  IF(O_i)
        IF_mat=matrix(nrow=p,ncol=N_trt)

        # [2022-12-09] George: We added this line to assing the 'u' object, which
        # we think is "user" that comes from the "elig_pid" column.
        u=sort(unique(baseline$elig_pid))
        for(i in 1:N_trt){
            # w=which(Visits_df_a$elig_pid==u[i])

            # if(length(w) ==0){ temp = 0 }
            # if(length(w) ==1){ temp=Term1_mat[,w] }
            # if(length(w) > 1){ temp= }
            IF_mat[,i]=rowSums(Term1_mat[,Visits_df_a$elig_pid==u[i], drop=FALSE]) + Term2_mat[,i]
        }


        ########  Target parameter (p x 1)
        Beta_hat_a=rowMeans(IF_mat)
        # means_fn_a=function(t, ...){
        #     B_t=matrix(spline_fn(t, ...),ncol=p)
        #     B_t%*%Beta_hat_a
        # }
        if(do.var){
            #! this can benefit from armadillo and matrix structure.
            Var_beta=
                (1/N_trt^2)*matrix(rowSums(apply(IF_mat - Beta_hat_a, 2, tcrossprod)),p,p)
            #(1/N_trt^2)*rowSums(Var_array,dims=2)
            curve_variance <- apply(B_t_matrix, 2, function(x)crossprod(x,  Var_beta %*% x))
        }
        purrr::compact(tibble(
            alpha = alpha,
            beta = list(Beta_hat_a),
            var.beta = list(Var_beta),
            mean_curve = list(tibble(
                time = spline_seq,
                mean = crossprod(B_t_matrix, Beta_hat_a)[,1],
                var = if(do.var) curve_variance
            ))
        ))
        # Means_a=matrix(sapply(spline_seq,means_fn_a, knots),ncol=1)
        # Means_mat[,(j+1)]=Means_a
    }

    V=tcrossprod(B_t_matrix)    #B_t_matrix%*%t(B_t_matrix)
    V_inverse=solve(V)

    Weights_term2=V_inverse%*%B_t_matrix

    p <- nrow(B_t_matrix)

    ###########################    Intensity model  ############################
    ######   Andersen-Gill model stratifying by assessment number
    # survival_trt = dplyr::filter(data$Survival, !!trt.expr)
    environment(intensity.formula) <- environment()
    AG_model = coxph(intensity.formula, id=survival[[patient_id_var]], data=survival)
    gamma = AG_model$coefficients


    ############  Estimated baseline intensities for each stratum  #############
    AG_surv=survfit(AG_model,newdata=data.frame(Prev_outcome=0))
    strata=AG_surv$strata


    base_intens <- dplyr::bind_cols(purrr::map2(AG_surv$strata, cumsum(dplyr::lag(AG_surv$strata, 1, 0)),
        function(stratum, cumprev){
            cumhaz <- data.frame( time   = AG_surv$time[(1:stratum)+cumprev]
                                , cumhaz = AG_surv$cumhaz[(1:stratum)+cumprev]
                                )
            sapply(1:End,lambda0_fn,b=30,surv=cumhaz) #< not sure if b is constant or should be passed as a parameter
        }))

    ###########################    Outcome model  ##############################

    Outcome_model <- glm.nb(outcome.formula,data = visits)


    theta = Outcome_model$theta

    coefs=Outcome_model$coefficients
    a=coefs[5]+coefs[6]

    Exp_gamma = exp(gamma * visits$Prev_outcome)
    mu_k = predict.glm(Outcome_model, newdata = visits, type='response')


    #########  Construct influence function for the spline parameter beta  #####
    Visits_df <- dplyr::filter( visits
                              , time >= min(spline_seq)
                              , time <= max(spline_seq)
                              ) |>
        dplyr::mutate(baseline_lambda = purrr::map2_dbl(time, Visit_number, ~base_intens[[.x, .y]]))

    K <- NROW(Visits_df)

    # baseline_lambda <- rep(NA,K)
    #
    # for(k in 1:K){
    #
    #     visit_number=Visits_df$Visit_number[k]
    #     time=Visits_df$time[k]
    #
    #     baseline_lambda[k]=base_intens[visit_number][time]
    # }


    purrr::map( Alpha_seq, estimate_curve_for_one_alpha, ...)
}
if(FALSE){#@Testing
    data <- ARC_data |>
        dplyr::filter(time <= 830) |>
        formatting_fn(
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
    spline_seq=seq(60,460,by=1)
    trt.expr <- rlang::expr(Trt == "home_visits")

    intensity.formula <- Surv(Prev_time,time,Event)~Prev_outcome+strata(Visit_number)
    outcome.formula   <- (6*Asthma_control)~ns(Prev_outcome,df=3) + time + Lag_time

    df <- estimate_mean_curves(
        data,
        treatment = Trt == "home_visits",
        patient_id_var = elig_pid,
        intensity.formula = Surv(Prev_time,time,Event)~Prev_outcome+strata(Visit_number),
        outcome.formula = (6*Asthma_control)~ns(Prev_outcome,df=3) + time + Lag_time,
        spline_seq = spline_seq,
        Knots = c( rep(min(spline_seq)-1, 4)
                 , mean(range(spline_seq))
                 , rep(max(spline_seq)+1, 4)
                 ),
        spline_fn = spline_fn,
        End = 830,
        Alpha_seq = c(-0.6,-0.3,0,0.3,0.6),
        Cond_mean_fn = Cond_mean_fn
    )

    dplyr::bind_rows(df) |>
        select(alpha, mean_curve) |>
        tidyr::unnest(mean_curve) |>
        ggplot2::qplot(data=_, ggplot2::aes(x=time, y=mean, col=alpha), geom='line')

}
if(F){
    knots <- c(59,59,59,59,260,461,461,461,461)
    spline_seq=seq(60,460,by=1)



    data <- ARC_data |>
        dplyr::filter(time <= 830) |>
        formatting_fn(
            id_var          = "elig_pid",
            treatment_var   = "Trt",
            outcome_var     = "Asthma_control",
            visitnumber_var = "Visit_number",
            time_var        = "time",
            V               = 5,
            knots           = knots,
            spline_seq      = spline_seq,
            last_day = 830
        )
    trt.expr <- rlang::expr(Trt == "home_visits")

    B_t_matrix <- matrix(sapply(spline_seq, spline_fn, knots),byrow=FALSE,nrow=length(knots)-4)



    mean_curve_ests(
        data$Baseline,
        data$Visits)

    result <- estimate_curve_for_one_alpha(
        alpha = 0
        , baseline = filter(data$Baseline, !!trt.expr)
        , survival = filter(data$Survival, !!trt.expr)
        , visits   = filter(data$Visits  , !!trt.expr)
        , patient_id_var = 'elig_pid'
        , B_t_matrix =B_t_matrix
        , spline_seq = spline_seq
        , intensity.formula =
            Surv(Prev_time,time,Event)~Prev_outcome+strata(Visit_number)
        , outcome.formula =
            (6*Asthma_control)~ns(Prev_outcome,df=3) + time + Lag_time
        , End = attr(data, 'End')
        , Cond_mean_fn = Cond_mean_fn
        , knots = attr(data, 'knots')
        , do.var=TRUE
    )

}
