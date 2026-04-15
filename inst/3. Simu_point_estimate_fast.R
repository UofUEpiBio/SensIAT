rm(list=ls())

library(dplyr)
library(survival)
library(MAVE)
library(dfoptim)
library(orthogonalsplinebasis) 
library(BB)

select <- dplyr::select


### Stage 1: - treatment group
{
  
  call_function_with_data <- function(index, fun, ...){
    
    data <- dataall[[index]][[1]] # for the original data
    fun(data, ...)
  }
  
  inner_ARC_SIR_analyze_data <- function(data, 
                                         N = 200, 
                                         End = 830,
                                         method = "dfsane"
  ){
    
    #  we use the seed from the external function 
    #  data is the dataset for analyzing: if we want to analyze the original data, we put df_sim
    #   if we want to analyze the bootstrapped data, we put bootstrapped data here
    df_sim <- data
    
    ####  period over which we are conducting inference
    spline_seq <- seq(60, 460, by=1)
    
    library(dplyr)
    library(survival)
    library(MAVE)
    library(dfoptim)
    library(orthogonalsplinebasis) 
    library(BB)
    
    select = dplyr::select
    
    
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
      #####   Single index model
      
      SIDR_Ravinew <- function(X, Y,
                               Y.CP = NULL,
                               initial = NULL,
                               kernel = "dnorm",
                               method = "optim",
                               optim_method = "BFGS", 
                               abs.tol = 1e-4,
                               bandwidth = NULL,
                               wi.boot = NULL,
                               index_ID)
      {
        X <- as.matrix(X)
        Y <- as.matrix(Y)
        
        number_n <- dim(X)[1]
        number_p <- dim(X)[2]
        
        if (is.null(initial))
        {
          initial <- c(1, rep(0, number_p-1))
        }else
        {
          initial <- as.vector(initial)
          initial <- initial/initial[1]
        }
        
        if (is.null(bandwidth))
        {
          if (kernel=="K2_Biweight")
          {
            if (is.null(wi.boot))
            {
              Eij3 <- function(parameter){
                K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h) 
                
                b <- c(1, parameter[1:(number_p-1)])
                h <- exp(parameter[number_p])
                
                x <- c(X%*%b) 
                y <- Y
                
                n <- length(y)
                yo <- order(y)
                ys <- y[yo]
                uy <- rle(ys)[[1]]
                cols <- cumsum(uy)  # Ravi 
                ei <- rep(0, n)
                for (i in 1:n){
                  Kih <- K(x-x[i],h=h)
                  
                  # remove all obs of ith obs's patient
                  index_remove <- which(index_ID == index_ID[i])
                  Kih[index_remove] <- 0
                  
                  # Kih[i] <- 0                       # the fix
                  denom <- sum(Kih)
                  ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                }
                return(sum(ei)/n^2)
              }
            }else
            {
              stop("There's no weighted version of the K2_Biweight kernel.")
            }
          }else if (kernel == "dnorm")
          {
            if (is.null(wi.boot))
            {
              Eij3 <- function(parameter){
                K <- function(x, h) dnorm(x/h, 0, 1) 
                
                b <- c(1, parameter[1:(number_p-1)])
                h <- exp(parameter[number_p])
                
                x <- c(X%*%b) 
                y <- Y
                
                n <- length(y)
                yo <- order(y)
                ys <- y[yo]
                uy <- rle(ys)[[1]]
                cols <- cumsum(uy)
                ei <- rep(0, n)
                for (i in 1:n){
                  Kih <- K(x-x[i],h=h)
                  
                  # remove all obs of ith obs's patient
                  index_remove <- which(index_ID == index_ID[i])
                  Kih[index_remove] <- 0
                  
                  # Kih[i] <- 0       # the fix
                  denom <- sum(Kih)
                  ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                }
                return(sum(ei)/n^2)
              }
            }else
            {
              stop("There's no weighted version of the dnorm kernel.")
            }
          }else if (kernel=="K4_Biweight")
          {
            if (is.null(wi.boot))
            {
              Eij3 <- function(parameter){
                K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h)
                
                b <- c(1, parameter[1:(number_p-1)])
                h <- exp(parameter[number_p])
                
                x <- c(X%*%b) 
                y <- Y
                
                n <- length(y)
                yo <- order(y)
                ys <- y[yo]
                uy <- rle(ys)[[1]]
                cols <- cumsum(uy)
                ei <- rep(0, n)
                for (i in 1:n){
                  Kih <- K(x-x[i],h=h)
                  
                  # remove all obs of ith obs's patient
                  index_remove <- which(index_ID == index_ID[i])
                  Kih[index_remove] <- 0
                  
                  # Kih[i] <- 0                       # the fix
                  denom <- sum(Kih)
                  ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                }
                return(sum(ei)/n^2)
              }
            }else
            {
              stop("There's no weighted version of the K4_Biweight kernel.")
            }
          }
          
          if(method == "nlminb")
          {
            esti <- nlminb(start = c(initial[-1], 0), 
                           objective = Eij3,
                           control = list(abs.tol = abs.tol))
          }else if (method == "optim")
          {
            # the new optimize function using optim, you can change the lower and upper
            esti <- optim(par = c(initial[-1], 0), 
                          fn = Eij3,
                          method = optim_method,
                          control = list(abstol = abs.tol))
          }else if (method == "nmk")
          {
            esti <- nmk(par = c(initial[-1], 0), 
                        fn = Eij3,
                        control = list(tol = abs.tol))
          }
          
          results <- list(coef = c(1, esti$par[1:(number_p-1)]),
                          bandwidth = exp(esti$par[number_p]),
                          details = esti)
        }else
        {
          if (kernel=="K2_Biweight")
          {
            if (is.null(wi.boot))
            {
              Eij3 <- function(parameter){
                K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h) 
                
                b <- c(1, parameter[1:(number_p-1)])
                h <- bandwidth
                
                x <- c(X%*%b) 
                y <- Y
                
                n <- length(y)
                yo <- order(y)
                ys <- y[yo]
                uy <- rle(ys)[[1]]
                cols <- cumsum(uy)  # Ravi 
                ei <- rep(0, n)
                for (i in 1:n){
                  Kih <- K(x-x[i],h=h)
                  
                  # remove all obs of ith obs's patient
                  index_remove <- which(index_ID == index_ID[i])
                  Kih[index_remove] <- 0
                  
                  # Kih[i] <- 0                       # the fix
                  denom <- sum(Kih)
                  ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                }
                return(sum(ei)/n^2)
              }
              
            }else
            {
              stop("There's no weighted version of the K2_Biweight kernel.")
            }
          }else if (kernel=="dnorm") 
          {
            if (is.null(wi.boot))
            {
              Eij3 <- function(parameter){
                K <- function(x, h) dnorm(x/h,0,1)
                
                b <- c(1, parameter[1:(number_p-1)])
                h <- bandwidth
                
                x <- c(X%*%b) 
                y <- Y
                
                n <- length(y)
                yo <- order(y)
                ys <- y[yo]
                uy <- rle(ys)[[1]]
                cols <- cumsum(uy)  # Ravi 
                ei <- rep(0, n)
                for (i in 1:n){
                  Kih <- K(x-x[i],h=h)
                  
                  # remove all obs of ith obs's patient
                  index_remove <- which(index_ID == index_ID[i])
                  Kih[index_remove] <- 0
                  
                  # Kih[i] <- 0                       # the fix
                  denom <- sum(Kih)
                  ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                }
                return(sum(ei)/n^2)
              }
            }else
            {
              stop("There's no weighted version of the dnorm kernel.")
            }
          }else if (kernel=="K4_Biweight")
          {
            if (is.null(wi.boot))
            {
              Eij3 <- function(parameter){
                K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h)
                
                b <- c(1, parameter[1:(number_p-1)])
                h <- bandwidth
                
                x <- c(X%*%b) 
                y <- Y
                
                n <- length(y)
                yo <- order(y)
                ys <- y[yo]
                uy <- rle(ys)[[1]]
                cols <- cumsum(uy)  # Ravi 
                ei <- rep(0, n)
                for (i in 1:n){
                  Kih <- K(x-x[i],h=h)
                  
                  # remove all obs of ith obs's patient
                  index_remove <- which(index_ID == index_ID[i])
                  Kih[index_remove] <- 0
                  
                  # Kih[i] <- 0                       # the fix
                  denom <- sum(Kih)
                  ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                }
                return(sum(ei)/n^2)
              }
            }else
            {
              stop("There's no weighted version of the K4_Biweight kernel.")
            }
          }
          
          if(method == "nlminb")
          {
            esti <- nlminb(start = initial[-1], 
                           objective = Eij3,
                           control = list(abs.tol = abs.tol))
          }else if (method == "optim")
          {
            # the new optimize function using optim, you can change the lower and upper
            esti <- optim(par = initial[-1], 
                          fn = Eij3,
                          method = optim_method,
                          control = list(abstol = abs.tol))
          }else if (method == "nmk")
          {
            esti <- nmk(par = initial[-1], 
                        fn = Eij3,
                        control = list(tol = abs.tol))
          }
          results <- list(coef = c(1, esti$par[1:(number_p-1)]),
                          bandwidth = bandwidth,
                          details = esti)
        }
        
        return(results)
      }
      
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
      
      # The general function in the package
      Cond_mean_fn_single <- function(alpha, X, Y, x, beta, bandwidth){
        
        y <- sort(unique(Y))
        # conditional distribution
        Fhat <- NW(Xb = X %*% beta, Y = Y,
                   xb = x %*% beta, y = y,
                   h = bandwidth)
        
        # density function
        Fhat1 <- c(0, Fhat[1:(length(y) - 1)]) 
        pmf <- Fhat - Fhat1
        
        E_exp_alphaY <- sum( exp(alpha*y)*pmf )
        
        E_Yexp_alphaY <- sum( y*exp(alpha*y)*pmf )
        
        E_Y_past <- E_Yexp_alphaY/E_exp_alphaY
        
        return(list(E_Y_past, E_exp_alphaY))
        
      }
      
      ##########################################
      
    }
    
    
    #############################################################
    #################   Weight Functions  #######################
    #############################################################
    {
      
      #####  the p x 401 matrix with columns B(t), t=60,...,460
      basis <- SplineBasis(knots = c(60,60,60,60,260,460,460,460,460))
      B_t_matrix <- evaluate(basis, spline_seq)
      
      # the length of beta
      p <- dim(B_t_matrix)[2]
      
      B3 <- evaluate(basis, 90)
      B6 <- evaluate(basis, 180)
      B9 <- evaluate(basis, 270)
      B12 <- evaluate(basis, 360)
      
      
      #####  we approximated the integral V=\int_t B(t)B(t)'dt using sums of rectangles of width 1
      V1 <- GramMatrix(basis)
      # Another methods to calculate V1
      # V1 <- t(B_t_matrix) %*% B_t_matrix 
      # V1 <- 0
      # for(i in 1:401){
      #   V1 <- V1 + as.matrix(t(B_t_matrix)[,i]) %*% B_t_matrix[i, ]
      # }
      V1_inverse <- solve(V1)
      
      # identity link
      Weights_term_identity <- V1_inverse%*%t(B_t_matrix)
      
      s_logit <- function(z){
        return(1/ (1 + exp(-z)))
      }
      # logit link - weight1
      Weights_term_logit1 <- function(t, beta){
        B_t <- evaluate(basis, t)
        return(V1_inverse %*% as.matrix(t(B_t) * as.vector(exp(B_t %*% beta)) * as.vector((1 + exp(-B_t %*% beta))^2)))
      }
      # logit link - weight2
      Weights_term_logit2 <- function(t, beta){
        B_t <- evaluate(basis, t)
        B_t_beta <- B_t_matrix %*% beta
        trans_B_t_beta <- exp(-B_t_beta) / (1 + exp(-B_t_beta))^2
        
        B_t_matrix_trans <- as.matrix(B_t_matrix * as.vector(trans_B_t_beta))
        V2 <- t(B_t_matrix) %*% B_t_matrix_trans
        # V2 <- 0
        # for(i in 1:nrow(B_t_matrix)){
        #   V2 <- V2 + as.matrix(t(B_t_matrix)[,i]) %*% B_t_matrix[i, ] * exp_B_t_beta[i]
        # }
        
        V2_inverse <- solve(V2)
        return(V2_inverse %*% t(B_t))
      } 
      
      
      s_log <- function(z){
        return(exp(z))
      }
      s_log_deriv <- function(z){
        return(exp(z))
      }
      # log link - weight1 
      Weights_term_log1 <- function(t, beta){
        B_t <- evaluate(basis, t)
        return(V1_inverse %*% as.matrix(t(B_t) * as.vector(exp(-B_t %*% beta))))
      }
      Weights_term_log1_deriv <- function(t, beta){
        B_t <- evaluate(basis, t)
        return(V1_inverse %*% as.matrix(t(B_t) * as.vector(exp(-B_t %*% beta))) %*% B_t )
      }
      
      
      # log link - weight2
      Weights_term_log2 <- function(t, beta){
        
        B_t <- evaluate(basis, t)
        B_t_beta <- B_t_matrix %*% beta
        exp_B_t_beta <- exp(B_t_beta)
        
        B_t_matrix_exp <- as.matrix(B_t_matrix * as.vector(exp_B_t_beta))
        V2 <- t(B_t_matrix) %*% B_t_matrix_exp
        # V2 <- 0
        # for(i in 1:nrow(B_t_matrix)){
        #   V2 <- V2 + as.matrix(t(B_t_matrix)[,i]) %*% B_t_matrix[i, ] * exp_B_t_beta[i]
        # }
        
        V2_inverse <- solve(V2)
        return(V2_inverse %*% t(B_t))
        
      }
    }
    
    
    #############################################################
    ############# Intensity and single index model  #############
    #############################################################
    {
      
      ###########  Fit intensity model to the simulated data
      {
        AG_model_sim <- coxph(Surv(prev_time,time,event) ~ Prev_outcome + strata(visit),
                              id = id, 
                              data = df_sim)
        # boot_id and id
        gamma <- AG_model_sim$coefficients
      }
      
      ############  Estimated baseline intensities
      {
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
      }
      
      ###########################    Outcome model - single index model  ###########################
      {
        
        All_visits_df_sim=filter(df_sim, event==1)%>%rename(Asthma_control=outcome)
        
        Xi <- data.frame(All_visits_df_sim$Prev_outcome,
                         All_visits_df_sim$time,
                         All_visits_df_sim$Lag_time)
        Xi <- as.matrix(Xi)
        Yi <- as.matrix(All_visits_df_sim$Asthma_control, ncol = 1)
        
        # new method to get the SIM's estimation
        initial <- coef(MAVE::mave.compute(Xi, Yi, max.dim = 1), 1)
        SID2 <- SIDR_Ravinew(X = Xi, 
                             Y = Yi,
                             initial = initial, 
                             kernel = "dnorm",
                             method = "nmk",
                             abs.tol = 1e-5,
                             index_ID = All_visits_df_sim$id)
        
      }
      
      ##################   Construct influence function for beta  ###############
      
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
      
    }
    
    
    #############################################################
    ################### Estimate beta  ##########################
    #############################################################
    # calculate the mean and variance for each alpha's value
    Alpha_seq <- c(-0.6, -0.3, 0, 0.3, 0.6)
    # Alpha_seq <- c(-0.6)
  
    
    #############################################################
    ################### Estimate beta's variance  ###############
    #############################################################
    Target_IF_var <- function(beta,
                              alpha){
      
      Visits_df_a <- Visits_df_sim
      
      #### Term1 of the influence function ####
      E_Y_past <- rep(NA,K)
      E_exp_alphaY <- rep(NA,K)
      Exp_gamma <- rep(NA,K)
      
      for(k in 1:K){
        
        df_k = Visits_df_a[k, ]
        Exp_gamma[k] = exp(gamma*df_k$Prev_outcome)
        
        temp <- Cond_mean_fn_single(alpha,
                                    X = Xi,
                                    Y = Yi,
                                    x = c(df_k$Prev_outcome,
                                          df_k$time,
                                          df_k$Lag_time),
                                    beta = SID2$coef,
                                    bandwidth=SID2$bandwidth)
        
        E_Y_past[k]=temp[[1]]
        E_exp_alphaY[k]=temp[[2]]
      }
      
      Visits_df_a <- mutate(Visits_df_a, Exp_gamma, E_Y_past, E_exp_alphaY)
      Visits_df_a <- mutate(Visits_df_a, 
                            Term1_unweighted = (Asthma_control - E_Y_past)/
                              (baseline_lambda*Exp_gamma*exp(-alpha*Asthma_control)*E_exp_alphaY) )
      
      ##########   the kth column is the term corresponding to visit k, in Term 1
      Term1_mat = matrix(nrow = p, ncol = K)
      Term1_mat_var = array(dim=c(p, p, K))
      
      for(k in 1:K){
        
        time_k <- Visits_df_a$time[k]
        Term1_mat[, k] =  Weights_term_log1(time_k, beta) * Visits_df_a$Term1_unweighted[k]
        Term1_mat_var[ , , k] = Weights_term_log1_deriv(time_k, beta) * Visits_df_a$Term1_unweighted[k]
        
      }
      
      
      ############  Term 2 of the influence function:
      ############   the ith column is Term 2 for participant i
      Term2_mat <- matrix(nrow = p, ncol = N)
      Term2_mat_var_1 = array(dim = c(p, p, N))
      Term2_mat_var_2 = array(dim = c(p, p, N))
      
      length_time <- length(spline_seq)
      
      # start <- Sys.time()
      for(i in 1:N){
        print(i)
        df_i1 <- filter(df_sim, id == i)
        
        ############## change this part for single index model ##############
        {
          # for the time interval spline_seq, generate the covariate matrix X
          df_est <- data.frame(time = spline_seq,
                               Lag_time = rep(0, length_time),
                               Prev_outcome =  rep(0, length_time))
          
          # new version
          {
            for(k in 1:length_time){
              t <- spline_seq[k]
              
              min_time <- min(df_i1$time)
              max_time <- max(df_i1$time)
              
              if(t %in% df_i1$time){
                index <- which(df_i1$time == t)
                df_est[k, 2:3] <- c(df_i1$Lag_time[index], 
                                    df_i1$Prev_outcome[index])
              }
              else if(t < min_time){
                df_est[k, 2:3] <- c(t, 
                                    df_i1$Prev_outcome[1])
              }
              else if(t > max_time){
                temp <- df_i1$outcome[nrow(df_i1)]
                df_est[k, 2:3] <- c(t - max_time,
                                    temp) 
              }
              else{
                t_index1 <- max(which(df_i1$time < t))
                temp <- df_i1$outcome[t_index1]
                df_est[k, 2:3] <- c(t - df_i1$time[t_index1], 
                                    temp)
              }
              
            }
          }  
          
          Time_means_single = matrix(0, nrow = p, ncol = length_time)
          Time_means_single_var1 = array(dim = c(p, p, length_time))
          Time_means_single_var2 = array(dim = c(p, p, length_time))
          
          for(l in 1:length_time){
            temp <- Cond_mean_fn_single(alpha,
                                        X = Xi,
                                        Y = Yi,
                                        x = c(df_est$Prev_outcome[l],
                                              df_est$time[l],
                                              df_est$Lag_time[l]),
                                        beta = SID2$coef,
                                        bandwidth = SID2$bandwidth)
            
            B_t_l <- evaluate(basis, df_est$time[l])
            s_t_l <- s_log(B_t_l %*% beta)
            W_t_l <- Weights_term_log1(df_est$time[l], beta)
            
            Time_means_single[, l] <- W_t_l * as.numeric((temp[[1]] - s_t_l))
            Time_means_single_var1[, , l] <- Weights_term_log1_deriv(df_est$time[l], beta) * as.numeric((temp[[1]] - s_t_l))
            Time_means_single_var2[, , l] <- W_t_l %*% as.matrix(as.numeric(s_log_deriv(B_t_l %*% beta)) * B_t_l)
            
          }
          
          # The following is the same formula as the package 
          Term2_mat[, i] <- (Time_means_single[, -401] %*% rep(1, length_time - 1) + 
                               Time_means_single[, -1] %*% rep(1, length_time - 1) )/2
          
          Term2_mat_var_1[, , i] <- (apply(Time_means_single_var1[, , -401], c(1, 2), sum) + 
                                       apply(Time_means_single_var1[, , -1], c(1, 2), sum))/2
          
          Term2_mat_var_2[, , i] <- (apply(Time_means_single_var2[, , -401], c(1, 2), sum) + 
                                       apply(Time_means_single_var2[, , -1], c(1, 2), sum))/2
          
        }
      }
      # end <- Sys.time()
      # end - start
      # 33 s
      
      #######   IF
      IF_mat <- matrix(nrow = p, ncol = N)
      Temp <- matrix(nrow = p, ncol = N)
      #######   Variance 
      Var_array_deriv <- array(dim = c(p, p, N))
      Temp_var <- array(dim = c(p, p, N))
      Var_array_phi  <- array(dim = c(p, p, N))

      for(i in 1:N){
        
        w = which(Visits_df_a$id == i)
        
        if(length(w) ==0){
          temp1 = 0
          temp2 = 0
        }
        if(length(w) ==1){
          temp1 = Term1_mat[, w]
          temp2 = Term1_mat_var[ , , w]
        }
        if(length(w) > 1){
          temp1 = rowSums(Term1_mat[, w])
          temp2 = apply(Term1_mat_var[ , , w], c(1, 2), sum)
        }
        
        Temp[, i] <- temp1
        IF_mat[,i] <- temp1 + Term2_mat[,i]
       
        Temp_var[ , , i] <- temp2
        Var_array_deriv[ , , i] <- Temp_var[ , , i] + Term2_mat_var_1[, , i] -  Term2_mat_var_2[, , i]
        
        temp <- IF_mat[, i]
        Var_array_phi[ , , i] <- temp%*%t(temp)
        
      }
      
      Var_deriv_inverse <- solve(rowMeans(Var_array_deriv, dims = 2))
      Var_beta_hat <- t(Var_deriv_inverse) %*% rowMeans(Var_array_phi, dims=2) %*% Var_deriv_inverse
      
      return(list(Beta = beta,
                  Var_beta = 1/N * Var_beta_hat)) 
      
    }
    
    Means_mat <- matrix(ncol=5, nrow=length(Alpha_seq))
    colnames(Means_mat) = c("Alpha", "month3", "month6","month9", "month12")
    Means_mat[, 1] <- Alpha_seq
    
    Vars_mat <- matrix(ncol=5, nrow=length(Alpha_seq))
    colnames(Vars_mat) = c("Alpha", "month3", "month6","month9", "month12")
    Vars_mat[, 1] <- Alpha_seq
    
    for(j in 1:length(Alpha_seq)){
      
      #### Prepare for term 1 and term 2 ####
      {
        #### Prepare for term 1 ####
        alpha <- Alpha_seq[j]
        
        Visits_df_a <- Visits_df_sim
        
        E_Y_past <- rep(NA,K)
        E_exp_alphaY <- rep(NA,K)
        Exp_gamma <- rep(NA,K)
        
        for(k in 1:K){
          
          df_k = Visits_df_a[k, ]
          Exp_gamma[k] = exp(gamma*df_k$Prev_outcome)
          
          temp <- Cond_mean_fn_single(alpha,
                                      X = Xi,
                                      Y = Yi,
                                      x = c(df_k$Prev_outcome,
                                            df_k$time,
                                            df_k$Lag_time),
                                      beta = SID2$coef,
                                      bandwidth=SID2$bandwidth)
          
          E_Y_past[k]=temp[[1]]
          E_exp_alphaY[k]=temp[[2]]
        }
        
        Visits_df_a <- mutate(Visits_df_a, Exp_gamma, E_Y_past, E_exp_alphaY)
        Visits_df_a <- mutate(Visits_df_a, 
                              Term1_unweighted = (Asthma_control - E_Y_past)/
                                (baseline_lambda*Exp_gamma*exp(-alpha*Asthma_control)*E_exp_alphaY) )
        
        #### Prepare for term 2 ####
        # Term2_mat <- matrix(nrow = p, ncol = N)
        length_time <- length(spline_seq)
        Time_means_single_all <- matrix(NA, nrow = length_time, ncol = N)
        
        # start <- Sys.time()
        for(i in 1:N){
          print(i)
          df_i1 <- filter(df_sim, id == i)
          
          ############## change this part for single index model ##############
          {
            # for the time interval spline_seq, generate the covariate matrix X
            df_est <- data.frame(time = spline_seq,
                                 Lag_time = rep(0, length_time),
                                 Prev_outcome =  rep(0, length_time))
            
            # new version
            {
              for(k in 1:length_time){
                t <- spline_seq[k]
                
                min_time <- min(df_i1$time)
                max_time <- max(df_i1$time)
                
                if(t %in% df_i1$time){
                  index <- which(df_i1$time == t)
                  df_est[k, 2:3] <- c(df_i1$Lag_time[index], 
                                      df_i1$Prev_outcome[index])
                }
                else if(t < min_time){
                  df_est[k, 2:3] <- c(t, 
                                      df_i1$Prev_outcome[1])
                }
                else if(t > max_time){
                  temp <- df_i1$outcome[nrow(df_i1)]
                  df_est[k, 2:3] <- c(t - max_time,
                                      temp) 
                }
                else{
                  t_index1 <- max(which(df_i1$time < t))
                  temp <- df_i1$outcome[t_index1]
                  df_est[k, 2:3] <- c(t - df_i1$time[t_index1], 
                                      temp)
                }
                
              }
            } 
            
            Time_means_single <- rep(0, length_time)
            # Time_means_single  <- matrix(NA, nrow = p, ncol = length_time)
            for(l in 1:length_time){
              temp <- Cond_mean_fn_single(alpha,
                                          X = Xi,
                                          Y = Yi,
                                          x = c(df_est$Prev_outcome[l],
                                                df_est$time[l],
                                                df_est$Lag_time[l]),
                                          beta = SID2$coef,
                                          bandwidth = SID2$bandwidth)
              
              Time_means_single[l] <- temp[[1]]
              # Time_means_single[, l] <- Weights_term_log1(df_est$time[l], beta) * as.numeric((temp[[1]] - s_log(evaluate(basis, df_est$time[l]) %*% beta)))
            }
            
            Time_means_single_all[, i] <-  Time_means_single
            # The following is the same formula as the package 
            # Term2_mat[, i] <- (Time_means_single[, -401] %*% rep(1, length_time - 1) + 
            #                      Time_means_single[, -1] %*% rep(1, length_time - 1) )/2
          }
        }
        # end <- Sys.time()
        # end - start
        # 33 s
      }
      
      Target_IF <- function(beta){
        
        ### Input parameter
        # Visits_df_a
        # Time_means_single_all
        # B_t_matrix
        # spline_seq
        # length_time
        # Weights_term_log1()
        # s_log()
        
        ##########   the kth column is the term corresponding to visit k, in Term 1
        Term1_mat = matrix(nrow=p, ncol=K)
        for(k in 1:K){
          time_k <- Visits_df_a$time[k]
          Term1_mat[,k] =  Weights_term_log1(time_k, beta) * Visits_df_a$Term1_unweighted[k]
        }
        
        ############  Term 2 of the influence function:
        ############   the ith column is Term 2 for participant i
        Term2_mat <- matrix(nrow = p, ncol = N)
        
        # start <- Sys.time()
        for(i in 1:N){
          # print(i)
          Time_means_single <- Time_means_single_all[ , i]
            
          ############## change this part for single index model ##############
          {
          
            Time_means_single_W  <- matrix(0, nrow = p, ncol = length_time)
            for(l in 1:length_time){
              Time_means_single_W[, l] <- Weights_term_log1(spline_seq[l], beta) * as.numeric((Time_means_single[l] - s_log(B_t_matrix[l,] %*% beta)))
            }
            
            # The following is the same formula as the package 
            Term2_mat[, i] <- (Time_means_single_W[, -401] %*% rep(1, length_time - 1) +
                                 Time_means_single_W[, -1] %*% rep(1, length_time - 1) )/2
          }
        }
        # end <- Sys.time()
        # end - start
        # 33 s

        #######   Subject-specific p x 1 influence function  IF(O_i) 
        IF_mat <- matrix(nrow=p, ncol=N)
        Temp <- matrix(nrow=p, ncol=N)
        
        for(i in 1:N){
          
          w = which(Visits_df_a$id == i)
          
          if(length(w) == 0){
            temp = 0
          }
          if(length(w) == 1){
            temp = Term1_mat[, w]
          }
          if(length(w) > 1){
            temp = rowSums(Term1_mat[, w])
          }
          
          Temp[,i] <- temp
          IF_mat[,i] = temp + Term2_mat[, i]
          
        }
        
        return(rowMeans(IF_mat)) 
        
      }
      
      if(method == "dfsane"){
        
        # set.seed(123)
        # start1 <- Sys.time()
        res1 <- dfsane(par = rep(0, p), fn = Target_IF, control = list(trace = FALSE))
        # end1 <- Sys.time()
        # end1 - start1
        
      }else if(method == "sane"){
        res1 <- sane(par = rep(0, p), fn = Target_IF, control = list(trace = FALSE))
      }else{
        print("No such method")
      }
      
    
      Est_res <- Target_IF_var(beta = res1$par,
                               alpha = Alpha_seq[j])
      
      ########  beta and the variance of beta
      Beta_hat <- Est_res$Beta
      Var_beta <- Est_res$Var_beta
      
      month3  <- s_log(B3 %*% Beta_hat)  
      Var_month3  <- (s_log(B3 %*% Beta_hat))^2 * B3 %*% Var_beta %*% t(B3)
      
      month6  <- s_log(B6 %*% Beta_hat) 
      Var_month6  <- (s_log(B6 %*% Beta_hat))^2 * B6 %*% Var_beta %*% t(B6)
      
      month9  <- s_log(B9 %*% Beta_hat)  
      Var_month9  <- (s_log(B9 %*% Beta_hat))^2 * B9 %*% Var_beta %*% t(B9)
      
      month12  <- s_log(B12%*%Beta_hat) 
      Var_month12  <- (s_log(B12 %*% Beta_hat))^2 * B12 %*% Var_beta %*% t(B12)
      
      Means_mat[j, 2:5] <- c(month3, month6, month9, month12)
      Vars_mat[j, 2:5] <- c(Var_month3, Var_month6, Var_month9, Var_month12)
      
      # Below is the previous identity link result
      
      # alpha = 0
      # month3   month6   month9  month12 
      # 1.863563 2.055067 1.958206 2.006252
      # month3      month6      month9     month12 
      # 0.009506188 0.007441910 0.007667211 0.008820933
      
      # alpha = -0.6
      # month3   month6   month9  month12 
      # 1.316642 1.511650 1.489761 1.499114
      # month3      month6      month9     month12 
      # 0.007706564 0.003268530 0.003994408 0.005164910
    }
    
    Results <- list(Means_mat, Vars_mat)
    
    return(Results) 
    
  }
  
}


setwd("/uufs/chpc.utah.edu/common/home/u6049227/R_output")
dataall <- readRDS('SIM_sim_data_treat')
# dataall <- readRDS('SIM_sim_data_control')

system.time(temp <- call_function_with_data(index = 1, fun = inner_ARC_SIR_analyze_data))
# 9.8 mins for all alpha values

library(parallel)

Index <- seq(401, 500, 1)
# Index <- 1

start <- Sys.time()
result_trt <- parallel::mclapply(Index, 
                                 call_function_with_data,
                                 inner_ARC_SIR_analyze_data,
                                 mc.cores = 64)
end <- Sys.time()
end - start

setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/ARC_simu_result")
saveRDS(result_trt, file = 'simu_treat/ARC_simu_res_trt_401_500_sane')
# saveRDS(result_trt, file = 'simu_control/ARC_simu_res_crl_401_500')
