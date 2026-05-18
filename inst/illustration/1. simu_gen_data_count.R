rm(list=ls())

library(dplyr)
library(survival)
library(MASS)
library(splines)
library(dfoptim)
library(orthogonalsplinebasis) # substitute R package (splines2) and function (bSpline)

select <- dplyr::select


call_function_with_seed <- function(seed, fun, ...){
  
  set.seed(seed)
  fun(...)
  
}


inner_ARC_SIR_sim_data_single <- function(N = 30,
                                          group = "control"){  
  
  # N=200       #####  size of simulated data for one arm
  
  ####  period over which we are conducting inference
  spline_seq <- seq(60, 460, by=1)
  
  
  library(dplyr)
  library(survival)
  library(MASS)
  library(splines)
  library(dfoptim)
  library(orthogonalsplinebasis) # substitute R package (splines2) and function (bSpline)
  
  select=dplyr::select
  
  
  #############################################################
  #################   Functions  ##############################
  #############################################################
  {
    
    ################################################################
    #####   This performs additional formatting on the cleaned ARC data
    #####   Input:  the cleaned ARC data & the day at which to place all censored observations
    #####   Output:  a list of 3 data frames:  baseline visits;  post-baseline visits;
    #####      and a data frame for use with the `survival' package, which has rows indicating 
    #####      censoring for subjects with fewer than 4 post-baseline assessments
    formatting_fn=function(df, last_day){
      
      K=dim(df)[1]
      
      Prev_outcome=rep(0,K)
      Prev_time=rep(0,K)
      
      for(i in 1:K){
        
        if(df$Visit_number[i]==0){Prev_outcome[i]=NA}
        if(df$Visit_number[i] !=0){Prev_outcome[i]=df$Asthma_control[i-1]}
        
      }
      
      for(i in 1:K){
        
        if(df$Visit_number[i]==0){Prev_time[i]=NA}
        if(df$Visit_number[i] !=0){Prev_time[i]=df$time[i-1]}
        
      }
      
      df=mutate(df,Prev_outcome,Prev_time)
      df=mutate(df,Lag_time=time-Prev_time,Event=1)
      
      
      ####  Each subject's baseline visit
      Baseline_df=filter(df,Visit_number==0)%>%
        select(-Event,-Prev_outcome,-Prev_time,-Lag_time)
      
      
      ####  Each subject's post-baseline visits
      Visits_df=filter(df,Visit_number !=0)%>%select(-Event)
      
      
      Survival_df=df
      
      u=sort(unique(Survival_df$elig_pid))
      H=length(u)
      
      for(h in 1:H){
        
        df_h=filter(Survival_df,elig_pid==u[h])
        v=dim(df_h)[1]
        
        if(v<5){
          
          temp_df=data.frame(elig_pid=u[h],Asthma_control=NA,
                             time=last_day,Trt=df_h$Trt[v],Visit_number=v,
                             Prev_outcome=df_h$Asthma_control[v],
                             Prev_time=df_h$time[v],Lag_time=last_day-df_h$time[v],
                             Event=0)
          
          Survival_df=rbind(Survival_df,temp_df)
          
        }
      }
      
      Survival_df=arrange(Survival_df,elig_pid)%>%filter(Visit_number != 0)
      
      return(list(Baseline_df,Visits_df,Survival_df))
      
    }
    ################################################
    
    
    ####################################################
    # Estimate the baseline intensity function in lambda(t,O(t)), the function if on page 12 
    
    #######    Function to kernel-smooth a cumulative intensity using an 
    #######        Epanechnikov kernel 
    #######    Inputs are time t, bandwidth b, and a 'survfit' object generated
    #######        by the 'survival' package
    lambda0_fn=function(t,b,surv){
      
      times=surv$time
      cumhaz=surv$cumhaz
      
      increment=rep(0,length(times))
      increment[1]=cumhaz[1]
      for(j in 2:length(times)){ increment[j]=cumhaz[j]-cumhaz[j-1] }
      
      f_t=(1/b)*sum( 0.75*(1-( (t-times)/b )^2 )*( abs(t-times)<b )*increment)
      return(f_t)
      
    }
    ####################################################
    
    
    ####################################################
    ######  Function to generate times of a given post-baseline assessment 
    ######  Uses Ogata's Thinning Algorithm to generate times that will have
    ######    a prescribed intensity -- in this case, the estimated intensity function
    ######    from the ARC data
    
    ######  Inputs are the Visit Number, vectors of subjects' previous assessment times 
    ######    & outcomes, a value lambda_star for the thinning, and 
    ######    baseline intensity vector lambda0_t and fitted cox model for desired 
    ######    intensity function
    ######  Outputs a data frame with the assessment time, event indicator, 
    ######     previous time & outcome, lag time, and visit number for each subject
    
    ######  Lambda_star should be a value > the maximum for the prescribed intensity
    Times_gen_fn=function(Start_vec=rep(0,N),
                          Visit=1,
                          Outcomes_vec, 
                          lambda0_t, 
                          lambda_star,
                          cox_model){ 
      
      Prev_outcome_df=data.frame(Prev_outcome=c(Outcomes_vec), 
                                 Visit_number = rep(Visit, length(Outcomes_vec)))
      
      predict_vec=predict(cox_model, newdata=Prev_outcome_df, type="lp",reference="zero")
      # "lp" means linear predictor
      # A reference of "zero" causes no centering to be done
      
      Times_vec=rep(NA,N)
      Event_vec=rep(1,N)
      
      for(i in 1:N){
        
        start=Start_vec[i]
        
        done=0
        while(done==0){
          
          ####  no kth visit if they had no (k-1)st visit
          if(Start_vec[i]==End+Visit-5){ 
            
            Times_vec[i]=End+Visit-4
            done=1
            
          }else{
            
            #####  generate a potential assessment time Tstar 
            S=rexp(1,rate=lambda_star)
            Tstar=ceiling(start+S)
            
            #####  if the candidate assessment time is too late, end with no additional
            #####     assessment time for that subject 
            if(Tstar >=  End+Visit-5){
              
              Times_vec[i]=End+Visit-4
              done=1
              
            }else{
              
              ######  if the candidate assessment time is not too late, accept it with 
              ######    a probability based on the the value of the intent
              lambda=lambda0_t[Tstar]*exp(predict_vec[i])
              
              U=runif(1)
              accept=1*(U < lambda/lambda_star)
              
              if(accept==1){
                
                Times_vec[i]=Tstar
                done=1
                
              }else{ start=Tstar }
              
            }
          } 
        }
        
        #####  Mark as having no new assessment time if applicable
        if( Times_vec[i]==(End+Visit-4) ){ Event_vec[i]=0 } 
        
      }
      
      df=data.frame(id=1:N,
                    time=Times_vec,
                    event=Event_vec,
                    Prev_outcome=Outcomes_vec,
                    prev_time=Start_vec,
                    Lag_time=Times_vec-Start_vec,
                    visit=Visit)
      
      return(df)
      
    }
    ###############################
    
    
    ###################################
    #####  Function for drawing a value of Y (support = 0, 1/6, 2/6, ..., 6)
    
    #####  Inputs:  a prediction from a negative binomial regression model on the 
    #####     'response' scale mu, and the theta (size) parameter from this model
    
    #####  Output:  1/6 times a draw from a negative binomial distribution with mean mu
    #####     and size theta, truncated to have support 0,1,2,..., 36
    NW <- function(Xb, Y, xb, y, h, kernel = "dnorm"){
      
      if(kernel == "dnorm"){
        K <- function(x, h){dnorm(x/h, 0, 1)} # Gaussian 
      } else if(kernel == "K2_Biweight"){
        K <- function(x, h){15/16*(1-(x/h)^2)^2 * (abs(x) <= h)} # K2_biweight
      } else if(kernel=="K4_Biweight"){
        K <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }# K4_biweight
      }
      
      Kxb <- sapply(xb, function(x, Xb) K(Xb-x, h), Xb=Xb)
      
      Ylty <- sapply(y, function(x, Y) 1*(Y <= x), Y=Y)
      
      denom <- colSums(Kxb)
      
      fyxb <- (denom!=0)*crossprod(Kxb, Ylty)/(denom + (denom==0))
      
      return(fyxb)
      
    }
    
    Y_draw_fn_single <- function(Fhat, y){
      
      # density function
      Fhat1 <- c(0, Fhat[1:(length(Fhat) - 1)]) 
      pmf <- Fhat - Fhat1
      
      pmf[which(pmf < 1e-10)] <- 0 
      pmf <- pmf / sum(pmf)
      
      Y_new <- sample(y, size = 1, replace = T, prob = pmf)
      return(Y_new)
      
    }
    
    ###################################
    
    f_transfer <- function(data){
      
      # ids <- unique(data$id)
      data2 <- data[0, ]
      
      for(i in 1:N){
        
        temp <- filter(data, id == i)
        temp1 <- filter(temp, !is.na(outcome))
        
        if(nrow(temp1) < 4 && nrow(temp1) >= 1){
          temp1 <- rbind(temp1,
                         c(i, End, 0, 
                           temp1$outcome[nrow(temp1)], # previous outcome
                           temp1$time[nrow(temp1)], # previous time
                           End - temp1$time[nrow(temp1)], # lag time
                           temp1$visit[nrow(temp1)] + 1, 
                           NA))
        }else if(nrow(temp1) == 0){
          print(i)
          temp1 <- temp[1, ]
          temp1$time <- End
          temp1$Lag_time <- End - temp$prev_time[1]
        }
        
        data2 <- rbind(data2, temp1)
        
      }
      
      return(data2)
      
    }
    
    
  }
  
  
  ############################################################
  ##############  Generate simulated data ##############
  ############################################################
  # simulate data with a sample size of n=200 for treatment group
  {
    ############################################################
    #######   Fit models to the ARC data
    ############################################################
    
    ARC_data <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_data.rds")
    
    y_unique <- unique(ARC_data$Asthma_control)
    y_unique <- y_unique[is.na(y_unique) == F] |> sort()
    
    End <- 830
    ARC_data2 = filter(ARC_data, time <= End)
    
    ARC_formatted=formatting_fn(df=ARC_data2,last_day = End)
    
    if(group == "control"){
      
      ARC_baseline=ARC_formatted[[1]]
      ARC_baseline_hv=filter(ARC_baseline, Trt=="control") # for control group
      n_hv=dim(ARC_baseline_hv)[1]
      
      ARC_visits=ARC_formatted[[2]]
      ARC_visits_hv=filter(ARC_visits, Trt=="control")
      
      ARC_survival=ARC_formatted[[3]]
      ARC_survival_hv=filter(ARC_survival, Trt=="control")
      
    }else{
      
      ARC_baseline=ARC_formatted[[1]]
      ARC_baseline_hv=filter(ARC_baseline,Trt=="home_visits") # for treatment group
      
      n_hv=dim(ARC_baseline_hv)[1]
      ARC_visits=ARC_formatted[[2]]
      ARC_visits_hv=filter(ARC_visits,Trt=="home_visits")
      
      ARC_survival=ARC_formatted[[3]]
      ARC_survival_hv=filter(ARC_survival,Trt=="home_visits")
      
    }
    
    ###########################    Intensity model  ######################
    
    AG_model_hv=coxph(Surv(Prev_time,time,Event) ~ Prev_outcome + strata(Visit_number),
                      id=elig_pid,
                      data=ARC_survival_hv)
    
    ############  Estimated baseline intensities ############
    {
      AG_surv_hv=survfit(AG_model_hv,newdata=data.frame(Prev_outcome=0))
      strata=AG_surv_hv$strata
      
      v1=strata[1]
      cumhaz_v1=data.frame(time=AG_surv_hv$time[1:v1],cumhaz=AG_surv_hv$cumhaz[1:v1])
      base_intens_v1=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)
      
      
      v2=strata[2]
      cumhaz_v2=data.frame(time=AG_surv_hv$time[(v1+1):(v1+v2)],
                           cumhaz=AG_surv_hv$cumhaz[(v1+1):(v1+v2)])
      base_intens_v2=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)
      
      
      v3=strata[3]
      cumhaz_v3=data.frame(time=AG_surv_hv$time[(v1+v2+1):(v1+v2+v3)],
                           cumhaz=AG_surv_hv$cumhaz[(v1+v2+1):(v1+v2+v3)])
      base_intens_v3=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)
      
      
      v4=strata[4]
      cumhaz_v4=data.frame(time=AG_surv_hv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
                           cumhaz=AG_surv_hv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
      base_intens_v4=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)
    }
    
    ######   Outcome model - single index model #######
    {
      
      # Outcome_model_hv <- glm.nb((6*Asthma_control)~Prev_outcome+time+Lag_time,
      #                            data=ARC_visits_hv)
      # 
      # theta_ARC=Outcome_model_hv$theta
      
      Xi <- data.frame(ARC_visits_hv$Prev_outcome,
                       ARC_visits_hv$time,
                       ARC_visits_hv$Lag_time)
      Xi <- as.matrix(Xi)
      Yi <- as.matrix(ARC_visits_hv$Asthma_control, ncol = 1)
      
      res2 <- SensIAT::fit_SensIAT_single_index_fixed_coef_model(formula = Asthma_control ~ -1 + Prev_outcome + time + Lag_time,
                                                                 data = ARC_visits_hv,
                                                                 kernel = "dnorm",
                                                                 method = "nmk",
                                                                 id = elig_pid, 
                                                                 abs.tol = 1e-7,
                                                                 initial = NULL)
      SIM_coef <- res2$coef
      SIM_band <- res2$bandwidth
      
      Xb = Xi %*% SIM_coef
      
    }
    
    #############################################################
    ##############     Generate simulated data    ###############
    #############################################################
    {
      
      samp = sample(n_hv, size = N, replace=TRUE)
      
      ################################################
      
      Y0_sim = ARC_baseline_hv[samp, ]
      Y0_sim = data.frame(Y0_sim)%>%select(Asthma_control)%>%
        rename(Prev_outcome=Asthma_control)
      
      df1_sim=Times_gen_fn(lambda0_t=base_intens_v1,
                           Outcomes_vec=Y0_sim,
                           lambda_star=0.02,
                           cox_model = AG_model_hv)
      
      newdat_1=data.frame(Prev_outcome=Y0_sim,
                          time=df1_sim$time,
                          Lag_time=df1_sim$Lag_time)
      
      ############## The single index model ##############
      {
        # generate the probability of Yi = 1 at the first assessment time
        Y1_sim <- rep(NA, N)
        for(n_index in 1:N){
          temp <- NW(Xb = Xb, Y = Yi,
                     xb = as.matrix(newdat_1[n_index, ]) %*% SIM_coef, y = y_unique,
                     h = SIM_band)
          Y1_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
        }
        
        df1_sim = mutate(df1_sim,outcome=Y1_sim)
        for(i in 1:N){
          if(df1_sim$event[i]==0){ df1_sim$outcome[i]=NA }
        }
      }
      
      
      ################################################
      
      df2_sim=Times_gen_fn(Start_vec=df1_sim$time, 
                           Visit=2,
                           lambda0_t=base_intens_v2,
                           Outcomes_vec=df1_sim$outcome,
                           lambda_star=0.02,
                           cox_model=AG_model_hv)
      
      newdat_2=data.frame(Prev_outcome=Y1_sim,
                          time=df2_sim$time,
                          Lag_time=df2_sim$Lag_time)
      
      ############## The single index model ##############
      {
        # generate the probability of Yi = 1 at the first assessment time
        Y2_sim <- rep(NA, N)
        for(n_index in 1:N){
          temp <- NW(Xb = Xb, Y = Yi,
                     xb = as.matrix(newdat_2[n_index, ]) %*% SIM_coef, y = y_unique,
                     h = SIM_band)
          Y2_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
        }
        
        df2_sim=mutate(df2_sim,outcome=Y2_sim)
        for(i in 1:N){
          if(df2_sim$event[i]==0){ df2_sim$outcome[i]=NA }
        }
      }
      
      
      ################################################
      
      df3_sim=Times_gen_fn(Start_vec=df2_sim$time, 
                           Visit=3,
                           lambda0_t=base_intens_v3,
                           Outcomes_vec=df2_sim$outcome,
                           lambda_star=0.02,
                           cox_model=AG_model_hv)
      
      newdat_3=data.frame(Prev_outcome=Y2_sim,
                          time=df3_sim$time,
                          Lag_time=df3_sim$Lag_time)
      
      ############## The single index model ##############
      {
        # generate the probability of Yi = 1 at the first assessment time
        Y3_sim <- rep(NA, N)
        for(n_index in 1:N){
          temp <- NW(Xb = Xb, Y = Yi,
                     xb = as.matrix(newdat_3[n_index, ]) %*% SIM_coef, y = y_unique,
                     h = SIM_band)
          Y3_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
        }
        
        df3_sim=mutate(df3_sim,outcome=Y3_sim)
        for(i in 1:N){
          if(df3_sim$event[i]==0){ df3_sim$outcome[i]=NA }
        }
      }
      
      
      ################################################
      
      df4_sim=Times_gen_fn(Start_vec=df3_sim$time, 
                           Visit=4,
                           lambda0_t=base_intens_v4,
                           Outcomes_vec=df3_sim$outcome,
                           lambda_star=0.02,
                           cox_model=AG_model_hv)
      
      newdat_4=data.frame(Prev_outcome=Y3_sim,
                          time=df4_sim$time,
                          Lag_time=df4_sim$Lag_time)
      
      ############## The single index model ##############
      {
        # generate the probability of Yi = 1 at the first assessment time
        Y4_sim <- rep(NA, N)
        for(n_index in 1:N){
          temp <- NW(Xb = Xb, Y = Yi,
                     xb = as.matrix(newdat_4[n_index, ]) %*% SIM_coef, y = y_unique,
                     h = SIM_band)
          Y4_sim[n_index] <- Y_draw_fn_single(Fhat = temp, y = y_unique)
        }
        
        df4_sim=mutate(df4_sim,outcome=Y4_sim)
        for(i in 1:N){
          if(df4_sim$event[i]==0){ df4_sim$outcome[i]=NA }
        }
      }
      
      ################################################
      
      df_sim_full <- rbind(df1_sim,df2_sim,df3_sim,df4_sim) %>% arrange(id,time)
      
    }
    
  }
  
  
  data_trans <- f_transfer(df_sim_full)
  
  
  for(i in 1:N){
    if(i == 1){
      temp <- filter(data_trans, id == i)
      # print(paste0(i,",", nrow(temp)))
      df_sim_temp <- data.frame(id = i, 
                                time = temp$prev_time[1], 
                                outcome = temp$Prev_outcome[1], 
                                Visit_number = 0)
      for(j in 1:nrow(temp)){
        df_sim_temp <- rbind(df_sim_temp,
                             data.frame(id = i, 
                                        time = temp$time[j], 
                                        outcome = temp$outcome[j], 
                                        Visit_number = temp$visit[j]))
      }
    }
    else{
      temp <- filter(data_trans, id == i)
      print(paste0(i,",", nrow(temp)))
      df_sim_temp1 <- data.frame(id = i, 
                                 time = temp$prev_time[1], 
                                 outcome = temp$Prev_outcome[1], 
                                 Visit_number = 0)
      for(j in 1:nrow(temp)){
        df_sim_temp1 <- rbind(df_sim_temp1,
                              data.frame(id = i, 
                                         time = temp$time[j], 
                                         outcome = temp$outcome[j], 
                                         Visit_number = temp$visit[j]))
      }
      
      df_sim_temp <- rbind(df_sim_temp, df_sim_temp1)
      
    }
  }
  
  df_sim_temp <- filter(df_sim_temp, is.na(outcome) == F)
  
  Results <- list(data_trans, df_sim_full, df_sim_temp)
  
  return(Results) 
  
}


sim_data <- inner_ARC_SIR_sim_data_single(N = 30, group = "treat")
saveRDS(sim_data, "/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/sim_data_30_trt.RData")




