
wald <- function(mean_orig, var_orig, Zstar = qnorm(0.975)){
  CI_left <- mean_orig - Zstar * sqrt(var_orig)
  CI_right <- mean_orig + Zstar * sqrt(var_orig)
  
  return(c(CI_left, CI_right))
}
logit <- function(x){
  return(log( x / (1-x)))
}
expit <- function(x){
  return(1 / (1 + exp(-x)))
}



### confidence interval for E[Y(t)]
paraboot_CI_EY <- function(Est_EY, # original mean estimate
                           Est_paraboot_EY, # parametric bootstrap mean estimates (\widehat{E[Y(t)]}_k)
                           Est_paraboot_Bt, # parametric bootstrap linear predictor estimates (\beta_k^\prime B(t))
                           link_type = "log" # log, logit and identity
                          ){

  if(link_type == "log"){
    
    result <- list(Wald(Est_EY, var(Est_paraboot_EY)),
                   Wald(log(Est_EY), var(Est_paraboot_Bt)) |> exp() )
  }
  
  if(link_type == "logit"){
    
    result <- list(Wald(Est_EY, var(Est_paraboot_EY)),
                   Wald(logit(Est_EY), var(Est_paraboot_Bt)) |> expit()  )
    
  }
  
  if(link_type == "identity"){
    result <- list(Wald(Est_EY, var(Est_paraboot_EY)))
  }
  
  
  return(result)
}


### confidence interval for the treatment effect
paraboot_CI_effect <- function(Est_EY_trt, # original mean estimate (treat)
                               Est_EY_crl, # original mean estimate (control)
                               Est_paraboot_EY_trt, # parametric bootstrap mean estimates (\widehat{E[Y(t)]}_k) (treat)
                               Est_paraboot_EY_crl, # parametric bootstrap mean estimates (\widehat{E[Y(t)]}_k) (control)
                               link_type = "log" # log, logit and identity
){
  
  Est_effect <- Est_EY_trt - Est_EY_crl
  Est_effect_EY <- Est_paraboot_EY_trt - Est_paraboot_EY_crl 
  result <- list(Wald(Est_effect, var(Est_effect_EY)))
  return(result)
  
}

