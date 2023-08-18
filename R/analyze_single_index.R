#' Analyze data with a single index model
#'
#' @param data
#' @param N Number of independent observational units, i.e. individuals, patients, etc., for single arm of simulated data.
#' @param spline_seq Evaluation points over which we are conducting inference
#' @param knots The knots for fitting the cublic spline.
#' @param time_mean
#' @param time_sd
#' @param Lag_time_mean
#' @param Lag_time_sd
#'
#' @return
#' @export
#'
#' @examples
analyze_single_index_model <-
function(
    data,
    N = 200,
    spline_seq=seq(60, 460, by=1),
    knots=c(59,59,59,59,260,461,461,461,461),
    End = NULL,

       time_mean = NULL,
       time_sd = NULL,
       Lag_time_mean = NULL,
       Lag_time_sd = NULL){

    # N=200       #####  size of simulated dataset for one arm
    # R=200       #####  number of bootstraps to be used

    #  we use the seed from the external function
    #  data is the dataset for analyzing: if we want to analyze the original data, we put df_sim
    #   if we want to analyze the bootstrapped data, we put bootstrapped data here
    df_sim <- data


    ####  knots for the cubic spline basis

    #####  length of target parameter
    p=length(knots)-4

    End=830

    # library(MYHRcpp)
    # library(parallel)
    # library(dfoptim)
    #
    # library(dplyr)
    # library(survival)
    # library(MASS)
    # select=dplyr::select




    ## Fit intensity model to the simulated data ----

            AG_model_sim <- coxph(Surv(prev_time,time,event)~Prev_outcome+strata(visit),
                                  id=id, data=df_sim)
            # boot_id and id
            gamma <- AG_model_sim$coefficients


    ## Estimated baseline intensities -----
            AG_surv_sim=survfit(AG_model_sim,newdata=data.frame(Prev_outcome=0))
            strata_sim=AG_surv_sim$strata


            v1_sim=strata_sim[1]
            cumhaz_v1_sim=data.frame(time=AG_surv_sim$time[1:v1_sim],
                                     cumhaz=AG_surv_sim$cumhaz[1:v1_sim])
            base_intens_v1_sim=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1_sim)


            v2_sim=strata_sim[2]
            cumhaz_v2_sim=data.frame(time=AG_surv_sim$time[(v1_sim+1):(v1_sim+v2_sim)],
                                     cumhaz=AG_surv_sim$cumhaz[(v1_sim+1):(v1_sim+v2_sim)])
            base_intens_v2_sim=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2_sim)


            v3_sim=strata_sim[3]
            cumhaz_v3_sim=data.frame(time=AG_surv_sim$time[(v1_sim+v2_sim+1):(v1_sim+v2_sim+v3_sim)],
                                     cumhaz=AG_surv_sim$cumhaz[(v1_sim+v2_sim+1):(v1_sim+v2_sim+v3_sim)])
            base_intens_v3_sim=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3_sim)


            v4_sim=strata_sim[4]
            cumhaz_v4_sim=data.frame(time=AG_surv_sim$time[(v1_sim+v2_sim+v3_sim+1):(v1_sim+v2_sim+v3_sim+v4_sim)],
                                     cumhaz=AG_surv_sim$cumhaz[(v1_sim+v2_sim+v3_sim+1):(v1_sim+v2_sim+v3_sim+v4_sim)])
            base_intens_v4_sim=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4_sim)

    ## Outcome model - single index model  ------
            All_visits_df_sim=filter(df_sim, event==1)%>%rename(Asthma_control=outcome)

            # Prev_outcome_mean <- mean(All_visits_df_sim$Prev_outcome)
            # Prev_outcome_sd <- sd(All_visits_df_sim$Prev_outcome)

            # ARC_visits_hv_SIM <- All_visits_df_sim
            Xi <- data.frame(All_visits_df_sim$Prev_outcome,
                             All_visits_df_sim$time_scale,
                             All_visits_df_sim$Lag_time_scale)
            Xi <- as.matrix(Xi)
            Yi <- as.matrix(All_visits_df_sim$Asthma_control, ncol = 1)

            # new method to get the SIM's estimation
            SDR1 <- cumuSIR(X = Xi, Y = Yi)
            SID2 <- SIDRnew(X = Xi, Y = Yi, initial = SDR1$basis[, 1], kernel = "dnorm", method = "nmk")

    ## Construct influence function for beta  -----

        ####  number of assessments in [a,b]
        Visits_df_sim <- filter(All_visits_df_sim, time >= min(spline_seq),time <= max(spline_seq))
        K <- dim(Visits_df_sim)[1]

        #######   Get each subject's baseline intensity at each of
        #######     their own visit times
        baseline_lambda=rep(NA,K)

        for(k in 1:K){

            visit_number=Visits_df_sim$visit[k]
            time=Visits_df_sim$time[k]

            if(visit_number==1){

                baseline_lambda[k]=base_intens_v1_sim[time]
            }

            if(visit_number==2){

                baseline_lambda[k]=base_intens_v2_sim[time]
            }

            if(visit_number==3){

                baseline_lambda[k]=base_intens_v3_sim[time]
            }

            if(visit_number==4){

                baseline_lambda[k]=base_intens_v4_sim[time]
            }

        }

        Visits_df_sim <- mutate(Visits_df_sim, baseline_lambda)

        # calculate the mean and variance for each alpha's value
        Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)

        Means_mat <- matrix(ncol=5, nrow=length(Alpha_seq))
        colnames(Means_mat) = c("Alpha", "month3", "month6","month9", "month12")
        Means_mat[, 1] <- Alpha_seq

        Vars_mat <- matrix(ncol=5, nrow=length(Alpha_seq))
        colnames(Vars_mat) = c("Alpha", "month3", "month6","month9", "month12")
        Vars_mat[, 1] <- Alpha_seq

        #####  the p x 441 matrix with columns B(t), t=60,...,460
        B_t_matrix <- matrix(sapply(spline_seq, spline_fn, knots=knots), byrow=FALSE, nrow=p)

        ####  approximating the integral using sums of rectangles of width 1
        V=B_t_matrix%*%t(B_t_matrix)
        V_inverse=solve(V)
        Weights_term2=V_inverse%*%B_t_matrix

        for(j in 1:length(Alpha_seq)){

            alpha <- Alpha_seq[j]
            print(alpha)

            E_Y_past <- rep(NA,K)
            E_exp_alphaY <- rep(NA,K)
            Exp_gamma <- rep(NA,K)

            # Xbeta <- Xi %*% SID2$coef

            ############## change this part for single index model ##############

            for(k in 1:K){

                df_k <- Visits_df_sim[k,]

                Exp_gamma[k] <- exp(gamma*df_k$Prev_outcome)

                temp <- Cond_mean_fn_single1(alpha,
                                             X = Xi,
                                             Y = Yi,
                                             x = c(df_k$Prev_outcome,
                                                   df_k$time_scale,
                                                   df_k$Lag_time_scale),
                                             beta=SID2$coef,
                                             bandwidth=SID2$bandwidth)

                E_Y_past[k]=temp[[1]]
                E_exp_alphaY[k]=temp[[2]]
            }

            Visits_df_sim_temp <- mutate(Visits_df_sim, Exp_gamma, E_Y_past, E_exp_alphaY)

            Visits_df_sim_temp <- mutate(Visits_df_sim_temp, Term1_unweighted = (Asthma_control-E_Y_past)/
                                             (baseline_lambda*Exp_gamma*exp(-alpha*Asthma_control)*E_exp_alphaY) )


            ##### Term 1 of the influence function
            #####  the kth column is the term corresponding to visit k
            Term1_mat=matrix(nrow=p,ncol=K)

            for(k in 1:K){

                time_k <- Visits_df_sim_temp$time[k]
                spline_k <- matrix(spline_fn(time_k, knots=knots), ncol=1)
                Term1_mat[,k] <- (V_inverse%*%spline_k)*Visits_df_sim_temp$Term1_unweighted[k]

            }


            #########  Term 2 of the influence function
            #########    the ith column is Term 2 for participant i
            Term2_mat = matrix(nrow=p,ncol=N)

            for(i in 1:N){
                print(i)
                df_i1 <- filter(df_sim, id==i)

                ############## change this part for single index model ##############
                {
                    # for the time interval spline_seq, generate the covariate matrix X
                    df_est <- data.frame(time = spline_seq,
                                         Lag_time = rep(0, length(spline_seq)),
                                         Prev_outcome =  rep(0, length(spline_seq)),
                                         time_scale = rep(0, length(spline_seq)),
                                         Lag_time_scale = rep(0, length(spline_seq)))

                    # new version
                    {
                        for(k in 1:length(spline_seq)){
                            t <- spline_seq[k]

                            min_time <- min(df_i1$time)
                            max_time <- max(df_i1$time)

                            if(t <= min_time){
                                df_est[k, 2:5] <- c(t,
                                                    df_i1$Prev_outcome[1],
                                                    (t - time_mean) / time_sd,
                                                    (t - Lag_time_mean)/Lag_time_sd)
                            }else if(t >= max_time){
                                temp <- df_i1$outcome[nrow(df_i1)]
                                df_est[k, 2:5] <- c(t - max_time,
                                                    temp,
                                                    (t - time_mean) / time_sd,
                                                    (t - max_time - Lag_time_mean) / Lag_time_sd)
                            }else{
                                # t within min_time and max_time
                                t_index1 <- max(which(df_i1$time < t))
                                t_index2 <- min(which(df_i1$time > t))
                                df_est[k, 2:5] <- c(t - df_i1$time[t_index1],
                                                    df_i1$Prev_outcome[t_index2],
                                                    (t - time_mean) / time_sd,
                                                    (t - df_i1$time[t_index1] - Lag_time_mean) / Lag_time_sd)
                            }

                        }
                    }

                    Time_means_single  <- rep(0, length(spline_seq))
                    for(l in 1:length(spline_seq)){
                        temp <- Cond_mean_fn_single1(alpha,
                                                     X = Xi,
                                                     Y = Yi,
                                                     x = c(df_est$Prev_outcome[l],
                                                           df_est$time_scale[l],
                                                           df_est$Lag_time_scale[l]),
                                                     beta=SID2$coef,
                                                     bandwidth=SID2$bandwidth)
                        Time_means_single[l] <- temp[[1]]
                    }
                    Term2_mat[,i] = Weights_term2%*%Time_means_single

                    ########## plug in something to break the loop if return value is NA ##########
                    if(T %in% is.na(Term2_mat[,i])){
                        print(i)
                        print("this bootstrap sample returns an error, we will run another bootstrap sample")
                    }
                }
            }

            #######   Subject-specific p x 1 influence function  IF(O_i)
            IF_mat=matrix(nrow=p,ncol=N)

            for(i in 1:N){

                w=which(Visits_df_sim$id==i)

                if(length(w) ==0){ temp = 0 }
                if(length(w) ==1){ temp=Term1_mat[,w] }
                if(length(w) > 1){ temp=rowSums(Term1_mat[,w]) }

                IF_mat[,i]=temp + Term2_mat[,i]

            }

            ########  Target parameter (p x 1)
            Beta_hat=rowMeans(IF_mat)

            ########  Influence function-based variance estimation
            Var_array=array(dim=c(p,p,N))

            for(i in 1:N){

                temp=IF_mat[,i]-Beta_hat

                Var_array[,,i]= temp%*%t(temp)

            }

            Var_beta_hat=(1/N^2)*rowSums(Var_array,dims=2)

            #######  Target-time means and IF-based variance estimates
            {
                B3=matrix(spline_fn(t=90,knots=knots),ncol=p)

                month_3=B3%*%Beta_hat

                Var_month3=B3%*%Var_beta_hat%*%t(B3)


                B6=matrix(spline_fn(t=180,knots=knots),ncol=p)

                month_6=B6%*%Beta_hat

                Var_month6=B6%*%Var_beta_hat%*%t(B6)


                B9=matrix(spline_fn(t=270,knots=knots),ncol=p)

                month_9=B9%*%Beta_hat

                Var_month9=B9%*%Var_beta_hat%*%t(B9)


                B12=matrix(spline_fn(t=360,knots=knots),ncol=p)

                month_12=B12%*%Beta_hat

                Var_month12=B12%*%Var_beta_hat%*%t(B12)

            }

            Means_mat[j, 2:5] <- c(month_3, month_6, month_9, month_12)
            Vars_mat[j, 2:5] <- c(Var_month3, Var_month6, Var_month9, Var_month12)

        }

    ## Output -----
        list(Means_mat, Vars_mat)
}
