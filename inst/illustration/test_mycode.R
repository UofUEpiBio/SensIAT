rm(list = ls())
setwd("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack")
# data <- readRDS('sim_data_30_trt.RData')
data <- readRDS('sim_data_30_trt_binary.RData')

df_sim <- data[[1]] # for the original data
colnames(df_sim)[4] <- "prev_outcome"
colnames(df_sim)[6] <- "delta_time"

df_intensity <- data.frame(id = df_sim$id,
                           time = df_sim$time,
                           outcome = df_sim$outcome,
                           visit = df_sim$visit + 1,
                           prev_outcome = df_sim$prev_outcome,
                           prev_time = df_sim$prev_time,
                           delta_time = df_sim$delta_time)
saveRDS(df_intensity, "data_intensity_mycode_binary")

start <- Sys.time()

N = 30
End = 830
method = "dfsane"
Alpha_seq = c(0)
SIM_model = "fixed_coef" # "norm1coef", "fixed_coef", "fixed_band"
SIM_tol = 1e-7
link_type = "logit" # "log", "logit", "identity"
IF_type = "weight1"
use_gauss_legendre = T
m_gl = 30 # Number of Gauss-Legendre nodes


library(dplyr)
library(survival)
library(MAVE)
library(dfoptim)
library(orthogonalsplinebasis) 
library(BB)
library(statmod)  # For Gauss-Legendre quadrature

select = dplyr::select

if (use_gauss_legendre) {
  # Gauss-Legendre nodes and weights
  gl <- gauss.quad(n = m_gl, kind = "legendre")
  # Transform from [-1, 1] to [60, 460]
  spline_seq <- 60 + (460 - 60) * (gl$nodes + 1) / 2
  w_gl <- gl$weights * (460 - 60) / 2  # Quadrature weights
  length_time <- m_gl
  cat("Using Gauss-Legendre quadrature with", m_gl, "nodes\n")
} else {
  # Original dense grid
  spline_seq <- seq(60, 460, by=1)
  length_time <- length(spline_seq)
  w_gl <- rep(1, length_time)  # For trapezoidal rule
  cat("Using original trapezoidal rule with", length_time, "points\n")
}

m <- length_time  

#  we use the seed from the external function 
#  data is the dataset for analyzing: if we want to analyze the original data, we put df_sim
#   if we want to analyze the bootstrapped data, we put bootstrapped data here
unique_ID <- unique(df_sim$id)


#############################################################
#################   Functions  ##############################
#############################################################
{
  
  
  ####################################################
  # Estimate the baseline intensity function in lambda(t,O(t)), the function if on page 12 
  
  #######    Function to kernel-smooth a cumulative intensity using an 
  #######        Epanechnikov kernel - OPTIMIZED VECTORIZED VERSION
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
  
  # Vectorized version of lambda0_fn
  lambda0_fn_vectorized = function(times, b, surv) {
    cumhaz = surv$cumhaz
    surv_times = surv$Time
    
    # Compute increments once
    increment = diff(c(0, cumhaz))
    
    # Create matrix for vectorized computation
    t_matrix = matrix(times, nrow = length(times), ncol = length(surv_times))
    surv_matrix = matrix(surv_times, nrow = length(times), ncol = length(surv_times), byrow = TRUE)
    inc_matrix = matrix(increment, nrow = length(times), ncol = length(surv_times), byrow = TRUE)
    
    # Epanechnikov kernel vectorized
    u = (t_matrix - surv_matrix) / b
    mask = abs(u) < 1
    kernel_vals = 0.75 * (1 - u^2) * mask
    
    # Sum over j dimension (columns)
    rowSums(kernel_vals * inc_matrix, na.rm = TRUE) / b
  }
  ####################################################
  
  
  ####################################################
  #####   Single index model
  
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
  
  # fast: handles matrix x (multiple points)
  Cond_mean_fn_single_optimized <- function(alpha, X, Y, x, beta, bandwidth){
    
    if (is.vector(x)) {
      x <- matrix(x, nrow = 1)
    }
    
    # Precompute Xb once
    Xb <- X %*% beta  # n_train x 1
    y_unique <- sort(unique(Y))
    
    Fhat <- NW(Xb = Xb, Y = Y,
               xb = x %*% beta,  # Handles matrix x
               y = y_unique,
               h = bandwidth)
    
    # Fhat is matrix: nrow(x) x length(y_unique)
    Fhat_prev <- cbind(0, Fhat[, -ncol(Fhat), drop = FALSE])
    pmf <- Fhat - Fhat_prev
    
    # Compute expectations with matrix operations
    y_mat <- matrix(y_unique, nrow = nrow(x), ncol = length(y_unique), byrow = TRUE)
    
    E_exp_alphaY <- rowSums(exp(alpha * y_mat) * pmf)
    E_Yexp_alphaY <- rowSums(y_mat * exp(alpha * y_mat) * pmf)
    E_Y_past <- E_Yexp_alphaY / E_exp_alphaY
    
    return(list(E_Y_past, E_exp_alphaY))
  }
  
  ##########################################
  
}


#############################################################
#################   Weight Functions  #######################
#############################################################
{
  
  #####  the p x m matrix with columns B(t) at quadrature nodes
  basis <- SplineBasis(knots = c(60,60,60,60,260,460,460,460,460))
  B_t_matrix <- evaluate(basis, spline_seq)
  
  p <- dim(B_t_matrix)[2]
  
  B3 <- evaluate(basis, 90)
  B6 <- evaluate(basis, 180)
  B9 <- evaluate(basis, 270)
  B12 <- evaluate(basis, 360)
  
  #####  we approximate the integral V=\int_t B(t)B(t)'dt using GL quadrature
  if (use_gauss_legendre) {
    V1 <- t(B_t_matrix) %*% diag(w_gl) %*% B_t_matrix
  } else {
    V1 <- GramMatrix(basis)
  }
  V1_inverse <- solve(V1)
  
  #############################################################
  #################   Vectorized Weight Functions  ############
  #############################################################
  
  # continuous Outcome; identity link
  {
    s_identity <- function(z) z
    s_identity_deriv <- function(z) 1
    
    compute_weights_vectorized_identity <- function(beta, type = "weight1") {
      return(V1_inverse %*% t(B_t_matrix))
    }
    
    compute_weights_deriv_vectorized_identity<- function(beta, weights, type = "weight1"){
      array(0, dim = c(p, p, m))
    }
    
    Weights_term_identity <- function(t, beta){
      # identity; weights1
      B_t <- evaluate(basis, t)
      return(V1_inverse %*% t(B_t))
    }
    
    Weights_term_identity_deriv <- function(t, beta){
      # identity; weights1; derivative
      return(matrix(0, nrow = p, ncol = p))
    }
    
  }
  
  # binary Outcome; logit link
  {
    
    s_logit <- function(x) 1 / (1 + exp(-x))
    s_logit_deriv <- function(x) exp(-x) / (1 + exp(-x))^2
    
    # Vectorized weight computation for ALL quadrature points at once
    compute_weights_vectorized_binary <- function(beta, type = "weight1") {
      B_beta <- as.vector(B_t_matrix %*% beta)
      
      if (type == "weight1") {
        ds_dz <- exp(B_beta) + 2 + exp(-B_beta)
        weights <- V1_inverse %*% t(B_t_matrix)  # p x m
        # Scale each column by ds_dz
        weights <- weights * matrix(ds_dz, nrow = p, ncol = m, byrow = TRUE)
      } else {
        ds_dz <- exp(B_beta) / (1 + exp(B_beta))^2
        V2 <- t(B_t_matrix) %*% diag(w_gl * ds_dz) %*% B_t_matrix
        V2_inv <- solve(V2)
        weights <- V2_inv %*% t(B_t_matrix)  # p x m
      }
      
      return(weights)
    }
    
    # Vectorized derivative computation
    compute_weights_deriv_vectorized_binary <- function(beta, weights, type = "weight1") {
      B_beta <- as.vector(B_t_matrix %*% beta)
      
      if (type == "weight1") {
        dds_dz <- exp(B_beta) - exp(-B_beta)
        derivs <- array(0, dim = c(p, p, m))
        
        for (l in 1:m) {
          B_l <- matrix(B_t_matrix[l, ], ncol = 1)  # p x 1
          
          # Compute outer product
          B_outer <- B_l %*% t(B_l)  # p x p
          
          # Derivative
          derivs[, , l] <- V1_inverse %*% B_outer * dds_dz[l]
        }
      } else {
        # For weight2 - more complex derivative
        ds_dz <- exp(B_beta) / (1 + exp(B_beta))^2
        dds_dz <- exp(B_beta) * (1 - exp(B_beta)) / (1 + exp(B_beta))^3
        
        # Compute V2
        V2 <- t(B_t_matrix) %*% diag(w_gl * ds_dz) %*% B_t_matrix
        V2_inv <- solve(V2)
        
        derivs <- array(0, dim = c(p, p, m))
        
        for (j in 1:m) {
          # scalar (B_k^T W2(t_j)) for all k: length m
          dot_k <- as.vector(B_t_matrix %*% weights[, j])
          
          # accumulate sum_k w_k * s''(eta_k) * (B_k^T W2(t_j)) * (B_k B_k^T)
          S <- matrix(0, p, p)
          for (k in 1:m) {
            Bk <- matrix(B_t_matrix[k, ], ncol = 1)  # p x 1
            S <- S + (w_gl[k] * dds_dz[k] * dot_k[k]) * (Bk %*% t(Bk))
          }
          
          derivs[, , j] <- - V2_inv %*% S
        }
        
      }
      return(derivs)
    }
    
    # logit link - weight1
    Weights_term_logit1 <- function(t, beta){
      # logit; weights1
      B_t <- evaluate(basis, t)
      return(V1_inverse %*% t(B_t) * as.vector(exp(B_t %*% beta) + 2 + exp(-B_t %*% beta)))
    }
    Weights_term_logit1_deriv <- function(t, beta){
      B_t <- evaluate(basis, t)  # 1 x p matrix
      B_t_vec <- as.vector(B_t)  # Convert to vector
      
      # Compute outer product
      B_outer <- outer(B_t_vec, B_t_vec)  # p x p
      
      return(V1_inverse %*% B_outer * as.vector(exp(B_t %*% beta) - exp(-B_t %*% beta)))
    }
    
    # logit link - weight2
    Weights_term_logit2 <- function(t, beta){
      B_t <- evaluate(basis, t)
      B_t_beta <- B_t_matrix %*% beta
      trans_B_t_beta <- as.vector(exp(B_t_beta) / (1 + exp(B_t_beta))^2)
      
      if (use_gauss_legendre) {
        # Gauss-Legendre weights
        B_t_matrix_trans <- B_t_matrix * trans_B_t_beta
        M1 <- t(B_t_matrix) %*% diag(w_gl * trans_B_t_beta) %*% B_t_matrix
      } else {
        # Original trapezoidal rule
        B_t_matrix_trans <- B_t_matrix * trans_B_t_beta
        M1 <- (t(B_t_matrix[-nrow(B_t_matrix), ]) %*% B_t_matrix_trans[-nrow(B_t_matrix), ] + 
                 t(B_t_matrix[-1, ]) %*% B_t_matrix_trans[-1, ]) / 2
      }
      
      M1_inverse <- solve(M1)
      
      return(M1_inverse %*% t(B_t))
    }
    Weights_term_logit2_deriv <- function(t, beta){
      B_t <- evaluate(basis, t)
      B_t_vec <- as.vector(B_t)
      B_t_beta <- B_t_matrix %*% beta
      trans_B_t_beta <- as.vector(exp(B_t_beta) / (1 + exp(B_t_beta))^2)
      trans2_B_t_beta <- as.vector(exp(B_t_beta) * (1 - exp(B_t_beta)) / (1 + exp(B_t_beta))^3)
      
      if (use_gauss_legendre) {
        V2 <- t(B_t_matrix) %*% diag(w_gl * trans_B_t_beta) %*% B_t_matrix
      } else {
        B_t_matrix_trans <- B_t_matrix * trans_B_t_beta
        V2 <- (t(B_t_matrix[-nrow(B_t_matrix), ]) %*% B_t_matrix_trans[-nrow(B_t_matrix), ] + 
                 t(B_t_matrix[-1, ]) %*% B_t_matrix_trans[-1, ]) / 2
      }
      
      V2_inv <- solve(V2)
      W2 <- V2_inv %*% t(B_t)  # p x 1
      B_outer <- outer(B_t_vec, B_t_vec)  # p x p
      
      if (use_gauss_legendre) {
        temp <- t(B_t_matrix) %*% diag(w_gl * trans2_B_t_beta * as.vector(B_t_matrix %*% W2)) %*% B_t_matrix
      } else {
        B_t_matrix_trans2 <- B_t_matrix * (trans2_B_t_beta * as.vector(B_t_matrix %*% W2))
        temp <- (t(B_t_matrix[-nrow(B_t_matrix), ]) %*% B_t_matrix_trans2[-nrow(B_t_matrix), ] +
                   t(B_t_matrix[-1, ]) %*% B_t_matrix_trans2[-1, ])/2
      }
      
      return(-V2_inv %*% temp)
    }
    
    
  }
  
  # count Outcome; log link
  {
    
    s_log <- function(z) exp(z)
    s_log_deriv <- function(z) exp(z)
    
    # Vectorized weight computation for ALL quadrature points at once
    compute_weights_vectorized_count <- function(beta, type = "weight1") {
      B_beta <- as.vector(B_t_matrix %*% beta)
      
      if (type == "weight1") {
        ds_dz <- exp(-B_beta) 
        weights <- V1_inverse %*% t(B_t_matrix)  # p x m
        # Scale each column by ds_dz
        weights <- weights * matrix(ds_dz, nrow = p, ncol = m, byrow = TRUE)
      } else {
        ds_dz <- exp(B_beta) 
        V2 <- t(B_t_matrix) %*% diag(w_gl * ds_dz) %*% B_t_matrix
        V2_inv <- solve(V2)
        weights <- V2_inv %*% t(B_t_matrix)  # p x m
      }
      
      return(weights)
    }
    
    # Vectorized derivative computation
    compute_weights_deriv_vectorized_count <- function(beta, weights, type = "weight1") {
      B_beta <- as.vector(B_t_matrix %*% beta)
      
      if (type == "weight1") {
        dds_dz <- - exp(-B_beta)
        derivs <- array(0, dim = c(p, p, m))
        
        for (l in 1:m) {
          B_l <- matrix(B_t_matrix[l, ], ncol = 1)  # p x 1
          
          # Compute outer product
          B_outer <- B_l %*% t(B_l)  # p x p
          
          # Derivative
          derivs[, , l] <- V1_inverse %*% B_outer * dds_dz[l]
        }
      } else {
        # For weight2 - more complex derivative
        ds_dz <- exp(B_beta) 
        dds_dz <- ds_dz
        
        # Compute V2
        V2 <- t(B_t_matrix) %*% diag(w_gl * ds_dz) %*% B_t_matrix
        V2_inv <- solve(V2)
        
        derivs <- array(0, dim = c(p, p, m))
        
        for (j in 1:m) {
          # scalar (B_k^T W2(t_j)) for all k: length m
          dot_k <- as.vector(B_t_matrix %*% weights[, j])
          
          # accumulate sum_k w_k * s''(eta_k) * (B_k^T W2(t_j)) * (B_k B_k^T)
          S <- matrix(0, p, p)
          for (k in 1:m) {
            Bk <- matrix(B_t_matrix[k, ], ncol = 1)  # p x 1
            S <- S + (w_gl[k] * dds_dz[k] * dot_k[k]) * (Bk %*% t(Bk))
          }
          
          derivs[, , j] <- - V2_inv %*% S
        }
        
      }
      return(derivs)
    }
    
    
    # log link - weight1 
    Weights_term_log1 <- function(t, beta){
      # log; weights1
      B_t <- evaluate(basis, t)
      return(V1_inverse %*% t(B_t) * as.vector(exp(-B_t %*% beta)))
    }
    Weights_term_log1_deriv <- function(t, beta){
      # log; weights1; derivative
      B_t <- evaluate(basis, t)
      return( - V1_inverse %*% t(B_t) %*% B_t * as.vector(exp(-B_t %*% beta)) )      
    }
    
    # log link - weight2
    Weights_term_log2 <- function(t, beta){
      B_t <- evaluate(basis, t)
      B_t_beta <- B_t_matrix %*% beta
      trans_B_t_beta <- as.vector(exp(B_t_beta))
      
      if (use_gauss_legendre) {
        # Gauss-Legendre weights
        B_t_matrix_trans <- B_t_matrix * trans_B_t_beta
        M1 <- t(B_t_matrix) %*% diag(w_gl * trans_B_t_beta) %*% B_t_matrix
      } else {
        # Original trapezoidal rule
        B_t_matrix_trans <- B_t_matrix * trans_B_t_beta
        M1 <- (t(B_t_matrix[-nrow(B_t_matrix), ]) %*% B_t_matrix_trans[-nrow(B_t_matrix), ] + 
                 t(B_t_matrix[-1, ]) %*% B_t_matrix_trans[-1, ]) / 2
      }
      
      M1_inverse <- solve(M1)
      
      return(M1_inverse %*% t(B_t))
    }
    Weights_term_log2_deriv <- function(t, beta){
      B_t <- evaluate(basis, t)
      B_t_vec <- as.vector(B_t)
      B_t_beta <- B_t_matrix %*% beta
      trans_B_t_beta <- as.vector(exp(B_t_beta))
      trans2_B_t_beta <- as.vector(exp(B_t_beta))
      
      if (use_gauss_legendre) {
        V2 <- t(B_t_matrix) %*% diag(w_gl * trans_B_t_beta) %*% B_t_matrix
      } else {
        B_t_matrix_trans <- B_t_matrix * trans_B_t_beta
        V2 <- (t(B_t_matrix[-nrow(B_t_matrix), ]) %*% B_t_matrix_trans[-nrow(B_t_matrix), ] + 
                 t(B_t_matrix[-1, ]) %*% B_t_matrix_trans[-1, ]) / 2
      }
      
      V2_inv <- solve(V2)
      W2 <- V2_inv %*% t(B_t)  # p x 1
      B_outer <- outer(B_t_vec, B_t_vec)  # p x p
      
      if (use_gauss_legendre) {
        temp <- t(B_t_matrix) %*% diag(w_gl * trans2_B_t_beta * as.vector(B_t_matrix %*% W2)) %*% B_t_matrix
      } else {
        B_t_matrix_trans2 <- B_t_matrix * (trans2_B_t_beta * as.vector(B_t_matrix %*% W2))
        temp <- (t(B_t_matrix[-nrow(B_t_matrix), ]) %*% B_t_matrix_trans2[-nrow(B_t_matrix), ] +
                   t(B_t_matrix[-1, ]) %*% B_t_matrix_trans2[-1, ])/2
      }
      
      return(-V2_inv %*% temp)
    }
    
  }
  
}


#############################################################
############# Intensity and single index model  #############
#############################################################
{
  
  # TIMING: Data preparation and model fitting
  cat("\n=== TIMING: Data preparation and model fitting ===\n")
  segment1_start <- Sys.time()
  
  ###########  Fit intensity model to the simulated data
  {
    AG_model_sim <- coxph(Surv(prev_time,time,event) ~ prev_outcome + strata(visit),
                          id = id, 
                          data = df_sim)
    # boot_id and id
    gamma <- AG_model_sim$coefficients
    
  }
  
  ############  Estimated baseline intensities - Vectorized
  {
    
    AG_surv_sim = survfit(AG_model_sim, newdata = data.frame(prev_outcome = 0))
    strata_sim = AG_surv_sim$strata
    
    # Create time sequence once
    all_times = 1:End
    
    # Split surv object by strata
    v1_sim = strata_sim[1]
    cumhaz_v1_sim = data.frame(Time = AG_surv_sim$time[1:v1_sim],
                               cumhaz = AG_surv_sim$cumhaz[1:v1_sim])
    
    v2_sim = strata_sim[2]
    cumhaz_v2_sim = data.frame(Time = AG_surv_sim$time[(v1_sim+1):(v1_sim+v2_sim)],
                               cumhaz = AG_surv_sim$cumhaz[(v1_sim+1):(v1_sim+v2_sim)])
    
    v3_sim = strata_sim[3]
    cumhaz_v3_sim = data.frame(Time = AG_surv_sim$time[(v1_sim+v2_sim+1):(v1_sim+v2_sim+v3_sim)],
                               cumhaz = AG_surv_sim$cumhaz[(v1_sim+v2_sim+1):(v1_sim+v2_sim+v3_sim)])
    
    v4_sim = strata_sim[4]
    cumhaz_v4_sim = data.frame(Time = AG_surv_sim$time[(v1_sim+v2_sim+v3_sim+1):(v1_sim+v2_sim+v3_sim+v4_sim)],
                               cumhaz = AG_surv_sim$cumhaz[(v1_sim+v2_sim+v3_sim+1):(v1_sim+v2_sim+v3_sim+v4_sim)])
    
    # Compute all at once - 4 calls instead of 3,320!
    base_intens_v1_sim = lambda0_fn_vectorized(all_times, b = 30, surv = cumhaz_v1_sim)
    base_intens_v2_sim = lambda0_fn_vectorized(all_times, b = 30, surv = cumhaz_v2_sim)
    base_intens_v3_sim = lambda0_fn_vectorized(all_times, b = 30, surv = cumhaz_v3_sim)
    base_intens_v4_sim = lambda0_fn_vectorized(all_times, b = 30, surv = cumhaz_v4_sim)
    
  }
  
  ###########################    Outcome model - single index model  ###########################
  {
    
    All_visits_df_sim = filter(df_sim, event==1)
    
    Xi <- data.frame(All_visits_df_sim$prev_outcome,
                     All_visits_df_sim$time,
                     All_visits_df_sim$delta_time)
    Xi <- as.matrix(Xi)
    Yi <- as.matrix(All_visits_df_sim$outcome, ncol = 1)
    
    if(SIM_model == "norm1coef"){
      res1 <- SensIAT::fit_SensIAT_single_index_norm1coef_model(formula = outcome ~ -1 + prev_outcome + time + delta_time,
                                                                data = All_visits_df_sim, 
                                                                kernel = "dnorm",
                                                                mave.method = "meanMAVE",
                                                                id = id,
                                                                bw.selection = "mse",
                                                                bw.method = "optim",
                                                                bw.range = c(0.01, 1.5),
                                                                reestimate.coef = 0)
      SIM_coef <- res1$coefficients
      SIM_band <- res1$bandwidth
      
    }else if(SIM_model ==  "fixed_coef"){
      
      res2 <- SensIAT::fit_SensIAT_single_index_fixed_coef_model(formula = outcome ~ -1 + prev_outcome + time + delta_time,
                                                                 data = All_visits_df_sim,
                                                                 kernel = "dnorm",
                                                                 method = "nmk",
                                                                 id = id, 
                                                                 abs.tol = SIM_tol,
                                                                 initial = NULL)
      SIM_coef <- res2$coef
      SIM_band <- res2$bandwidth
      
    }else{
      
      res3 <- SensIAT::fit_SensIAT_single_index_fixed_bandwidth_model(formula = outcome ~ -1 + prev_outcome + time + delta_time,
                                                                      data = All_visits_df_sim,
                                                                      kernel = "dnorm",
                                                                      method = "nmk",
                                                                      id = id, 
                                                                      abs.tol = SIM_tol,
                                                                      initial = NULL)
      SIM_coef <- res3$coefficients
      SIM_band <- res3$bandwidth
      
    }
    
  }
  
  ##################   Construct influence function for beta  ###############
  
  ####  number of assessments in [a,b]
  Visits_df_sim <- filter(All_visits_df_sim, time >= min(spline_seq), time <= max(spline_seq))
  K <- dim(Visits_df_sim)[1]
  
  #######   Get each subject's baseline intensity at each of 
  #######     their own visit times
  baseline_lambda = rep(NA, K)
  for(k in 1:K){
    
    visit_number = Visits_df_sim$visit[k]
    Time = Visits_df_sim$time[k]
    
    if(visit_number == 1){
      baseline_lambda[k] = base_intens_v1_sim[Time]
    } else if(visit_number == 2){
      baseline_lambda[k] = base_intens_v2_sim[Time]
    } else if(visit_number == 3){
      baseline_lambda[k] = base_intens_v3_sim[Time]
    } else if(visit_number == 4){
      baseline_lambda[k] = base_intens_v4_sim[Time]
    }
    
  }
  Visits_df_sim <- mutate(Visits_df_sim, baseline_lambda)
  
  segment1_end <- Sys.time()
  timing_results <- as.numeric(difftime(segment1_end, segment1_start, units = "secs"))
  cat(sprintf("Data preparation completed in %.2f seconds\n", timing_results))
  
}


#############################################################
################### Estimate beta  ##########################
#############################################################
# calculate the mean and variance for each alpha's value
length_alpha <- length(Alpha_seq)

#############################################################
################### Estimate beta's variance  ###############
#############################################################

Means_mat <- matrix(ncol=5, nrow=length(Alpha_seq))
colnames(Means_mat) = c("Alpha", "month3", "month6","month9", "month12")
Means_mat[, 1] <- Alpha_seq

Vars_mat <- matrix(ncol=5, nrow=length(Alpha_seq))
colnames(Vars_mat) = c("Alpha", "month3", "month6","month9", "month12")
Vars_mat[, 1] <- Alpha_seq


#### Prepare for term 1 and term 2 for all alpha ####
{
  
  Visits_df_a <- Visits_df_sim
  
  #### Term1 of the influence function ####
  E_Y_past_all <- matrix(NA, nrow = K, ncol = length_alpha)
  E_exp_alphaY_all <- matrix(NA, nrow = K, ncol = length_alpha)
  
  Exp_gamma <- rep(NA, K)
  
  Visits_df_a_all <- list()
  
  cat("Collecting all visit data for batch processing...\n")
  x_k_matrix <- matrix(NA, nrow = K, ncol = 3)
  for(k in 1:K) {
    df_k = Visits_df_a[k, ]
    Exp_gamma[k] = exp(gamma * df_k$prev_outcome)
    x_k_matrix[k, ] <- c(df_k$prev_outcome, df_k$time, df_k$delta_time)
  }
  
  # Batch compute for each alpha - 5 calls instead of K x 5
  cat("Batch computing conditional means for all visits...\n")
  for(j in 1:length_alpha) {
    alpha <- Alpha_seq[j]
    batch_results <- Cond_mean_fn_single_optimized(alpha,
                                                   X = Xi, Y = Yi,
                                                   x = x_k_matrix,  # All K visits
                                                   beta = SIM_coef,
                                                   bandwidth = SIM_band)
    
    E_Y_past_all[, j] <- batch_results[[1]]
    E_exp_alphaY_all[, j] <- batch_results[[2]]
  }
  
  for(j in 1:length_alpha){
    
    alpha <- Alpha_seq[j]
    
    Visits_df_a_j <- mutate(Visits_df_a, 
                            Exp_gamma, 
                            E_Y_past = E_Y_past_all[, j], 
                            E_exp_alphaY = E_exp_alphaY_all[, j])
    
    Visits_df_a_j <- mutate(Visits_df_a_j, 
                            Term1_unweighted = (outcome - E_Y_past) /
                              (baseline_lambda * Exp_gamma * exp(-alpha * outcome) * E_exp_alphaY) ) 
    
    Visits_df_a_all <- c(Visits_df_a_all, list(Visits_df_a_j))
    
  }
  
  #### Prepare for term 2
  Time_means_single_all <- lapply(1:length_alpha, function(i) matrix(NA, nrow = m, ncol = N))
  
  W_length <- rep(0, N)
  id_list <- list()
  
  all_x_matrix <- matrix(NA, nrow = N * m, ncol = 3)
  index_mapping <- matrix(NA, nrow = N * m, ncol = 2)  # (i, l)
  
  counter <- 1
  for(i in 1:N) {
    if (i %% 10 == 0) cat("Processing individual", i, "of", N, "\n")
    
    # For the influence function
    w <- which(Visits_df_a$id == unique_ID[i])
    W_length[i] <- length(w)
    id_list <- c(id_list, list(w))
    
    ############## Generate covariate matrix for ALL quadrature points for this individual ##############
    {
      df_i1 <- filter(df_sim, id == unique_ID[i])
      df_est <- data.frame(time = spline_seq,
                           delta_time = rep(0, m),
                           prev_outcome = rep(0, m))
      
      min_time <- min(df_i1$time)
      max_time <- max(df_i1$time)
      df_i1_outcome_max <- df_i1$outcome[nrow(df_i1)]
      df_i1_prev_outcome_1 <- df_i1$prev_outcome[1]
      
      for(l in 1:m) {
        t <- spline_seq[l]
        
        if(t %in% df_i1$time) {
          index <- which(df_i1$time == t)
          df_est[l, 2:3] <- c(df_i1$delta_time[index], df_i1$prev_outcome[index])
        } else if(t < min_time) {
          df_est[l, 2:3] <- c(t, df_i1_prev_outcome_1)
        } else if(t <= max_time) {
          t_index1 <- max(which(df_i1$time < t))
          temp <- df_i1$outcome[t_index1]
          df_est[l, 2:3] <- c(t - df_i1$time[t_index1], temp)
        } else {
          temp <- df_i1_outcome_max
          df_est[l, 2:3] <- c(t - max_time, temp) 
        }
      }
      
      for(l in 1:m) {
        all_x_matrix[counter, ] <- c(df_est$prev_outcome[l], df_est$time[l], df_est$delta_time[l])
        index_mapping[counter, ] <- c(i, l)
        counter <- counter + 1
      }
    }
  }
  
  # compute All conditional means in batches
  cat("Computing conditional means in batches...\n")
  
  for(j in 1:length_alpha) {
    alpha <- Alpha_seq[j]
    
    # Compute for ALL x at once - 1 call per alpha instead of N*m calls
    batch_results <- Cond_mean_fn_single_optimized(alpha, X = Xi, Y = Yi,
                                                   x = all_x_matrix,  # All quadrature points
                                                   beta = SIM_coef,
                                                   bandwidth = SIM_band)
    
    E_Y_past_batch <- batch_results[[1]]  # Length = N*m
    
    # Fill Time_means_single_all[[j]] using index_mapping
    for(idx in 1:(N*m)) {
      i <- index_mapping[idx, 1]
      l <- index_mapping[idx, 2]
      Time_means_single_all[[j]][l, i] <- E_Y_past_batch[idx]
    }
  }
  
}


#### Optimization process 
{
  
  # Helper function for trapezoidal integration
  compute_integral_trap <- function(weights_matrix, diff_vec){
    m <- length(diff_vec)
    temp <- weights_matrix[, -m] %*% rep(1, m - 1) + weights_matrix[, -1] %*% rep(1, m - 1)
    return(temp / 2)
  }
  
  if(link_type == "identity"){
    
    s_link <- s_identity
    s_link_deriv <- s_identity_deriv
    
    Weights_term = Weights_term_identity 
    Weights_term_deriv = Weights_term_identity_deriv
    Weights_term_vectorized = function(beta) compute_weights_vectorized_identity(beta, "weight1")
    Weights_term_deriv_vectorized = function(beta, weights) compute_weights_deriv_vectorized_identity(beta, weights, "weight1")
    
  }else if(link_type == "logit"){
    
    s_link <- s_logit
    s_link_deriv <- s_logit_deriv
    
    if(IF_type == "weight1"){
      Weights_term = Weights_term_logit1 
      Weights_term_deriv = Weights_term_logit1_deriv
      Weights_term_vectorized = function(beta) compute_weights_vectorized_binary(beta, "weight1")
      Weights_term_deriv_vectorized = function(beta, weights) compute_weights_deriv_vectorized_binary(beta, weights, "weight1")
    }else if(IF_type == "weight2"){
      Weights_term = Weights_term_logit2
      Weights_term_deriv = Weights_term_logit2_deriv
      Weights_term_vectorized = function(beta) compute_weights_vectorized_binary(beta, "weight2")
      Weights_term_deriv_vectorized = function(beta, weights) compute_weights_deriv_vectorized_binary(beta, weights, "weight2")
    }
    
  }else if(link_type == "log"){
    
    s_link <- s_log
    s_link_deriv <- s_log_deriv
    
    if(IF_type == "weight1"){
      Weights_term = Weights_term_log1 
      Weights_term_deriv = Weights_term_log1_deriv
      Weights_term_vectorized = function(beta) compute_weights_vectorized_count(beta, "weight1")
      Weights_term_deriv_vectorized = function(beta, weights) compute_weights_deriv_vectorized_count(beta, weights, "weight1")
    }else if(IF_type == "weight2"){
      Weights_term = Weights_term_log2
      Weights_term_deriv = Weights_term_log2_deriv
      Weights_term_vectorized = function(beta) compute_weights_vectorized_count(beta, "weight2")
      Weights_term_deriv_vectorized = function(beta, weights) compute_weights_deriv_vectorized_count(beta, weights, "weight2")
    }
    
  }
  
  # Optimized root-finding function
  Target_IF_fast_optimized <- function(beta){
    
    ### Input parameter
    # Visits_df_a_j
    # Time_means_single_all_j
    # B_t_matrix
    # spline_seq
    # length_time
    # Weights_term
    # s_link
    
    ##########   the kth column is the term corresponding to visit k, in Term 1
    Term1_mat = matrix(nrow = p, ncol = K)
    for(k in 1:K){
      
      time_k <- Visits_df_a_j$time[k]
      Term1_mat[,k] = Weights_term(time_k, beta) * Visits_df_a_j$Term1_unweighted[k]
      
    }
    
    ############  Term 2 of the influence function:
    Term2_mat <- matrix(0, nrow = p, ncol = N)
    
    # vectorized weight computation
    Weights_spline_seq <- Weights_term_vectorized(beta)  # p x m matrix
    
    # Precompute s(B(t)_beta) at all quadrature nodes
    B_beta <- B_t_matrix %*% beta  # m x 1
    s_B_t_beta <- as.vector(s_link(B_beta))  # length m
    
    # vectorized operations
    for(i in 1:N){
      
      Time_means_single <- Time_means_single_all_j[, i]
      diff_vec <- Time_means_single - s_B_t_beta  # length m
      
      if (use_gauss_legendre) {
        Term2_mat[, i] <- Weights_spline_seq %*% (w_gl * diff_vec)  # p x 1
      } else {
        # trapezoidal rule
        Time_means_single_W <- Weights_spline_seq * matrix(diff_vec, nrow = p, ncol = m, byrow = TRUE)
        Term2_mat[, i] <- compute_integral_trap(Time_means_single_W, diff_vec)
      }
    }
    
    #######   Subject-specific p x 1 influence function  IF(O_i) 
    IF_mat <- matrix(nrow = p, ncol = N)
    for(i in 1:N){
      
      w <- id_list[[i]]
      w_length <- W_length[i]
      
      if(w_length > 1){
        temp = rowSums(Term1_mat[, w])
      }else if(w_length == 1){
        temp = Term1_mat[, w]
      }else if(w_length == 0){
        temp = 0
      }
      
      IF_mat[, i] = temp + Term2_mat[, i]
      
    }
    
    return(rowMeans(IF_mat)) 
  }
  
  # Vectorized variance function
  Target_IF_var_fast_vectorized <- function(beta, alpha) {
    
    ### Precompute weights and derivatives
    weights <- Weights_term_vectorized(beta)  # p x m
    weight_derivs <- Weights_term_deriv_vectorized(beta, weights)  # p x p x m
    
    ### Compute common quantities
    B_beta <- B_t_matrix %*% beta  # m x 1
    s_vec <- s_link(B_beta)  # m x 1
    s_deriv_vec <- s_link_deriv(B_beta)  # m x 1
    
    ### Compute diff matrix for ALL individuals at once
    diff_matrix <- Time_means_single_all_j - matrix(s_vec, nrow = m, ncol = N)
    
    ### Compute Term2 for ALL individuals (vectorized)
    w_diff <- w_gl * diff_matrix  # m x N
    Term2_mat <- weights %*% w_diff  # p x N
    
    ### Compute Term2 derivatives (optimized)
    Term2_deriv1 <- array(0, dim = c(p, p, N))
    Term2_deriv2 <- array(0, dim = c(p, p, N))
    
    w_s_deriv <- w_gl * s_deriv_vec  # m x 1
    
    for (l in 1:m) {
      W_l <- weights[, l, drop = FALSE]  # p x 1
      B_l <- B_t_matrix[l, , drop = FALSE]  # 1 x p
      dW_l <- weight_derivs[, , l]  # p x p
      
      scalar_vector <- w_gl[l] * diff_matrix[l, ]  # length N
      
      temp_deriv1 <- array(0, dim = c(p, p, N))
      for (i in 1:N) {
        temp_deriv1[, , i] <- dW_l * scalar_vector[i]
      }
      
      Term2_deriv1 <- Term2_deriv1 + temp_deriv1
      outer_W_B <- W_l %*% B_l  # p x p
      Term2_deriv2_common <- w_s_deriv[l] * outer_W_B
      Term2_deriv2 <- Term2_deriv2 + array(Term2_deriv2_common, dim = c(p, p, N))
    }
    
    ### Term1 calculations
    Term1_mat <- matrix(nrow = p, ncol = K)
    Term1_deriv <- array(0, dim = c(p, p, K))
    
    for(k in 1:K) {
      time_k <- Visits_df_a_j$time[k]
      term1_k <- Visits_df_a_j$Term1_unweighted[k]
      
      W_k <- Weights_term(time_k, beta)
      dW_k <- Weights_term_deriv(time_k, beta)
      
      Term1_mat[, k] <- W_k * term1_k
      Term1_deriv[, , k] <- dW_k * term1_k
    }
    
    ### Combine terms for each individual
    IF_mat <- matrix(0, nrow = p, ncol = N)
    Term1 <- matrix(0, nrow = p, ncol = N)
    IF_deriv_total <- matrix(0, nrow = p, ncol = p)
    IF_square <- matrix(0, nrow = p, ncol = p)
    
    for(i in 1:N) {
      w <- id_list[[i]]
      w_length <- W_length[i]
      
      if(w_length > 0) {
        if(w_length > 1) {
          term1_sum <- rowSums(Term1_mat[, w, drop = FALSE])
          term1_deriv_sum <- apply(Term1_deriv[, , w, drop = FALSE], c(1, 2), sum)
        } else {
          term1_sum <- Term1_mat[, w]
          term1_deriv_sum <- Term1_deriv[, , w]
        }
      } else {
        term1_sum <- 0
        term1_deriv_sum <- 0
      }
      
      Term1[, i] <- term1_sum
      IF_mat[, i] <- term1_sum + Term2_mat[, i]
      
      IF_deriv_total <- IF_deriv_total + 
        term1_deriv_sum / N + 
        Term2_deriv1[, , i] / N - 
        Term2_deriv2[, , i] / N
      
      IF_square <- IF_square + tcrossprod(IF_mat[, i]) / N
    }
    
    ### Compute variance
    IF_deriv_inverse <- solve(IF_deriv_total)
    Var_beta_hat <- IF_deriv_inverse %*% IF_square %*% t(IF_deriv_inverse)
    
    return(list(
      Beta = beta,
      Var_beta = (1 / N) * Var_beta_hat,
      Term1, 
      Term2_mat
      # IF_deriv_inverse = IF_deriv_inverse,
      # IF = IF_mat
    ))
  }
  
  #############################################################
  ################### Main Analysis Loop  #####################
  #############################################################
  
  for(j in 1:length(Alpha_seq)){
    
    alpha <- Alpha_seq[j]
    Visits_df_a_j <- Visits_df_a_all[[j]]
    Time_means_single_all_j <- Time_means_single_all[[j]]
    
    cat("\n=== Solving for beta with alpha =", alpha, "===\n")
    
    # Root-finding
    rootfinding_start <- Sys.time()
    if(method == "dfsane"){
      res1_fast <- dfsane(par = rep(0, p), fn = Target_IF_fast_optimized, 
                          control = list(trace = FALSE, maxit = 1000))
    }else if(method == "sane"){
      res1_fast <- sane(par = rep(0, p), fn = Target_IF_fast_optimized, 
                        control = list(trace = FALSE, maxit = 1000))
    }else if(method == "BBsolve"){
      res1_fast <- BBsolve(par = rep(0, p), fn = Target_IF_fast_optimized)
    }else{
      stop("No such method")
    }
    cat("  Convergence:", res1_fast$convergence, "\n")
    
    # Variance estimation
    Est_res_fast <- Target_IF_var_fast_vectorized(beta = res1_fast$par, alpha = alpha)
    
    ########  beta and the variance of beta
    Beta_hat <- Est_res_fast$Beta
    Var_beta_hat <- Est_res_fast$Var_beta
    Term1 <-  Est_res_fast[[3]]
    Term2 <- Est_res_fast[[4]]
  }
  
}   

Results <- list(Beta_hat, Term1, Term2)
end <- Sys.time()

end - start
# Time difference of 3.780271 secs

save.image("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/mycode_output_logit.RData")
