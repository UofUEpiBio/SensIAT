###################################
#####  Function for drawing a value of Y (support = 0, 1/6, 2/6, ..., 6)

#####  Inputs:  a prediction from a negative binomial regression model on the
#####     'response' scale mu, and the theta (size) parameter from this model

#####  Output:  1/6 times a draw from a negative binomial distribution with mean mu
#####     and size theta, truncated to have support 0,1,2,..., 36
Y_draw_fn=function(mu,theta){

    pmf_fn=function(y){

        f_y=dnbinom(y,size=theta,prob=theta/(theta+mu))
        return(f_y)

    }
    pmf=matrix(sapply(0:36,pmf_fn),ncol=1)

    y=seq(0,6,by=1/6)
    Y=sample(y,size=1,prob=pmf)

    return(Y)

}
