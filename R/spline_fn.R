


#' Compute Spline Basis
#'
#' Computes the cubic spline basis \eqn{B(t)=(B_1(t), \dots, B_p(t))}
#' for a given day, `t`, using a given set of knots.
#'
#' @param t time
#' @param knots spline knots.
#'
#' @return vector of spline basis functions evaluated at `t`.
#'
#' @export
#'
#' @examples
#' spline_fn(260, c(59,59,59,59,260,461,461,461,461))
spline_fn <- function(t, knots=Knots){

    L=length(knots)

    #####  Iteratively construct the degree 3 spline basis, B3, from B0, B1, B2
    B0=rep(0,(L-1))
    for(i in 1:(L-1)){

        B0[i]=1*(t >=knots[i] & t< knots[i+1])

    }

    #################
    B1=rep(0,(L-2))

    B1[1]=0
    B1[2]=0
    B1[3]=(knots[5]-t)/(knots[5]-knots[4])*B0[4]
    B1[(L-4)]=(t-knots[(L-4)])/(knots[(L-3)]-knots[(L-4)])*B0[(L-4)]
    B1[(L-3)]=0
    B1[(L-2)]=0

    for(i in 4:(L-5)){

        B1[i]=(t-knots[i])/(knots[i+1]-knots[i])*B0[i]+(knots[i+2]-t)/(knots[i+2]-knots[i+1])*B0[i+1]

    }

    #################
    B2=rep(0,(L-3))

    B2[1]=0
    B2[2]=(knots[2+3]-t)/(knots[2+3]-knots[2+1])*B1[3]
    B2[(L-4)]=(t-knots[(L-4)])/(knots[(L-2)]-knots[(L-4)])*B1[(L-4)]
    B2[(L-3)]=0

    for(i in 3:(L-5)){

        B2[i]=(t-knots[i])/(knots[i+2]-knots[i])*B1[i]+(knots[i+3]-t)/(knots[i+3]-knots[i+1])*B1[i+1]

    }

    ###################
    B3=rep(0,(L-4))

    B3[1]=(knots[1+4]-t)/(knots[1+4]-knots[1+1])*B2[2]
    B3[(L-4)]=(t-knots[(L-4)])/(knots[(L-1)]-knots[(L-4)])*B2[(L-4)]

    for(i in 2:(L-5)){

        B3[i]=(t-knots[i])/(knots[i+3]-knots[i])*B2[i]+(knots[i+4]-t)/(knots[i+4]-knots[i+1])*B2[i+1]

    }

    return(B3)

}
