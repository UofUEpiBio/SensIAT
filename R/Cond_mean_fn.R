#' Estimate Subject Predicted Mean Outcome
#'
#' This function is used to get a subject's predicted mean outcome under a given
#' sensitivity parameter value `alpha`.
#' It is specific to the outcome model that we used on the
#' ARC data (which was negative binomial), and specific to the
#' tilting function alpha*Y(t), so would need to be changed if
#' using a different outcome model or different tilting function.
#'
#' @param mu mean
#' @param theta size
#' @param alpha sensitivity
#'
#' @return
#' a list containing
#' \eqn{\frac{E[ Y exp(\alpha Y) ]}{E[ exp(\alpha Y)]}}{(E[ Y*exp(\alpha Y) ]/(E[exp(\alpha Y)])}
#' and $E[ exp(\alpha Y) ]$
#' for Y a (truncated version of) a negative binomial having mean `mu` and
#' size `theta`.
#'
#' @export
#'
#' @examples
Cond_mean_fn <- function(mu,theta,alpha){

    pmf_fn=function(z){

        f_z=dnbinom(z,size=theta,prob=theta/(theta+mu))
        return(f_z)

    }

    #####  pmf truncated at z=36
    pmf=matrix(sapply(0:36,pmf_fn),ncol=1)

    const=sum(pmf)
    pmf=pmf/const

    y=seq(0,6,by=1/6)

    E_exp_alphaY=sum( exp(alpha*y)*pmf )

    E_Yexp_alphaY=sum( y*exp(alpha*y)*pmf )

    E_Y_past=E_Yexp_alphaY/E_exp_alphaY

    return(list(E_Y_past,E_exp_alphaY))
}
