rm(list = ls())
library(dplyr)
library(SensIAT)
library(survival)

setwd("/uufs/chpc.utah.edu/common/home/u6049227/SensIAT")
devtools::load_all()

data_all <- readRDS('/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/sim_data_30_trt.RData')
# data_all <- readRDS('/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/sim_data_30_trt_binary.RData')

sim_data <- data_all[[3]] # this is for R package, data_all[[1]] is for my code  
sim_data$Visit_number <- sim_data$Visit_number + 1
colnames(sim_data) <- c("id", "time", "outcome", "visit")

data_with_lags <- sim_data %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(
    prev_outcome = dplyr::lag(outcome, default = NA_real_, order_by = time),
    prev_time  = dplyr::lag(time, default = 0, order_by = time),
    delta_time  = time - dplyr::lag(time, default = NA_real_, order_by = time)
  ) %>%
  dplyr::ungroup()

# Fit models
# data_intensity_pack <- data_with_lags %>% dplyr::filter(time > 0) 
# This is different from my code. There is no NA outcome in the data_intensity_pack
data_intensity <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/data_intensity_mycode")
# data_intensity <- readRDS("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/data_intensity_mycode_binary")


intensity.model <- survival::coxph(
  Surv(prev_time, time, !is.na(outcome)) ~  prev_outcome  + strata(visit),
  id = id,
  data = data_intensity
)
attr(intensity.model, 'bandwidth') <- 30

# intensity.model


outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
  # Outcome ~ splines::ns( prev_outcome , df = 2) +  delta_time  - 1,
  outcome ~  prev_outcome  + time +  delta_time  - 1,
  id = id,
  data = data_with_lags %>% dplyr::filter(time > 0),
  kernel = "dnorm",
  abs.tol = 1e-7
)


impute_fn <- function(t, df) {
  data_wl <- df %>%
    dplyr::mutate(
       prev_time  = time,
       prev_outcome  = outcome,
       delta_time  = 0
    )
  extrapolate_from_last_observation(
    t, data_wl, "time",
    slopes = c("delta_time" = 1)
  )
}

knots <- c(60, 260, 460)


# Fit with different methods
system.time(
result_fast <- fit_SensIAT_marginal_mean_model_generalized(
  data = data_with_lags,
  time = data_with_lags$time,
  id = data_with_lags$id,
  alpha = 0,
  knots = knots,
  outcome.model = outcome.model,
  intensity.args = list(bandwidth = 30, 
                        kernel = "K2_weight"),
  intensity.model = intensity.model,
  loss = "lp_mse",
  link = "log",
  impute_data = impute_fn,
  term2_method = "fast"
)
)
# user  system elapsed 
# 0.498   0.000   0.987 



# check term1's calculation before the V1_inverse
{
  
  N <- 30
  Term1 <- matrix(NA, nrow = 5, ncol = N)
  for(i in 1:N){
    Term1[, i] <- result_fast$influence[[1]]$term1[[i]]
  }
  V1_inverse <- result_fast$V.inverse
  V1 <- solve(V1_inverse)
  Term1_new <- V1 %*% Term1 
  
}


mycode = new.env()
load("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/mycode_output.RData", envir = mycode)
# load("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/mycode_output_log.RData", envir = mycode)
# load("/uufs/chpc.utah.edu/common/home/u6049227/General_sensitivity/test_R_pack/mycode_output_logit.RData", envir = mycode)
mycode$Beta_hat
result_fast$coefficients


mycode$Term1 - Term1
mycode$AG_model_sim
intensity.model
all.equal(mycode$res2, outcome.model)

# mycode$res2$coef - outcome.model$coef
# mycode$res2$bandwidth - outcome.model$bandwidth


system.time(
result_original <- fit_SensIAT_marginal_mean_model_generalized(
  data = data_with_lags,
  time = data_with_lags$time,
  id = data_with_lags$id,
  alpha = 0,
  knots = knots,
  outcome.model = outcome.model,
  intensity.args = list(bandwidth = 30, 
                        kernel = "K2_weight"),
  intensity.model = intensity.model,
  loss = "lp_mse",
  link = "identity",
  impute_data = impute_fn,
  term2_method = "original",
  # term2_method = c("fast", "original", "fixed_grid", "seeded_adaptive", "gauss_legendre"),
  tol = 1e-12,  # Tight tolerance for convergence
  debug = TRUE  # Add debug flag if supported
  )
)
# user  system elapsed 
# 0.483   0.000   0.614
result_original$coefficients
mycode$Beta_hat


system.time(
result_fixed_grid <- fit_SensIAT_marginal_mean_model_generalized(
  data = data_with_lags,
  time = data_with_lags$time,
  id = data_with_lags$id,
  alpha = 0,
  knots = knots,
  outcome.model = outcome.model,
  intensity.args = list(bandwidth = 30, 
                        kernel = "K2_weight"),
  intensity.model = intensity.model,
  loss = "lp_mse",
  link = "identity",
  impute_data = impute_fn,
  term2_method = "fixed_grid",
  term2_grid_n = 200  # High density for accuracy
 )
)
# user  system elapsed 
# 12.121   0.012  12.734
result_fixed_grid$coefficients
mycode$Beta_hat


# check term1's calculation before the V1_inverse
{
  Term1 <- matrix(NA, nrow = 5, ncol = 10)
  for(i in 1:10){
    Term1[, i] <- result_fixed_grid$influence[[1]]$term1[[i]]
  }
  V1_inverse <- result_fixed_grid$V.inverse
  V1 <- solve(V1_inverse)
  Term1_new <- V1 %*% Term1 
}


system.time(
result_seeded <- fit_SensIAT_marginal_mean_model_generalized(
  data = data_with_lags,
  time = data_with_lags$time,
  id = data_with_lags$id,
  alpha = 0,
  knots = knots,
  outcome.model = outcome.model,
  intensity.args = list(bandwidth = 30, 
                        kernel = "K2_weight"),
  intensity.model = intensity.model,
  loss = "lp_mse",
  link = "identity",
  impute_data = impute_fn,
  term2_method = "seeded_adaptive",
  term2_grid_n = 100
)
)
# user  system elapsed 
# 40.685   0.041  41.658 
result_seeded$coefficients
mycode$Beta_hat


system.time(
  result_gauss_legendre <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$time,
    id = data_with_lags$id,
    alpha = 0,
    knots = knots,
    outcome.model = outcome.model,
    intensity.args = list(bandwidth = 30, 
                          kernel = "K2_weight"),
    intensity.model = intensity.model,
    loss = "lp_mse",
    link = "identity",
    impute_data = impute_fn,
    term2_method = "gauss_legendre",
    # term2_method = c("fast", "original", "fixed_grid", "seeded_adaptive", "gauss_legendre"),
    term2_grid_n = 200  # High density for accuracy
  )
)
# user  system elapsed 
# 18.305   0.120  19.227 
result_gauss_legendre$coefficients
mycode$Beta_hat

