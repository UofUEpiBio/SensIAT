#' Construct the Baseline Intensity
#'
#' @param visits.data  The Visits data formatted as the output from [formatting_fn]
#'      `Visits` component.
#' @param spline_seq The spline sequence.
#'
#' @return The `visists.data` filtered to those time points that fall within the
#'      spline_seq range with an added column `baseline_lambda` the subject's
#'      baseline intensity at each of their own visit times.
#' @export
#'
#' @examples
#' example("formatting_fn", package='pcoriRPackage')
#'
#' Visits_df <- ARC_formatted$Visits |>
#'     filter(, Trt == 'home_visits') |>
#'     construct_influence_function(60:460)
#' glimpse(Visits_df)
construct_baseline_intensity <-
function(visits.data, spline_seq){
    Visits_df <- filter(visits.data, time >= min(spline_seq),time <= max(spline_seq))
    K=dim(Visits_df)[1]

    baseline_lambda=rep(NA,K)

    for(k in 1:K){

        visit_number=Visits_df$Visit_number[k]
        time=Visits_df$time[k]


        if(visit_number==1){

            baseline_lambda[k]=base_intens_v1[time]
        }

        if(visit_number==2){

            baseline_lambda[k]=base_intens_v2[time]
        }

        if(visit_number==3){

            baseline_lambda[k]=base_intens_v3[time]
        }

        if(visit_number==4){

            baseline_lambda[k]=base_intens_v4[time]
        }

    }

    mutate(Visits_df, baseline_lambda)
}


#' Construct Influence Function
#'
#' For a given sensitivity parameter $\alpha$, first estimate $\beta$,
#' then estimate $Y(t)$.
#'
#' In the output `Term1_unweighted` is the estimator for
#' \deqn{
#'   \frac{
#'      Y_i(t_k)-\hat{E}[Y(t_k)|\overline{O}(t_k)_i]
#'   }{
#'      \hat{\rho}(t_k,Y_i(t_k),\overline{O}(t_k)_i)
#'   }.
#' }
#'
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
construct_influence_function <-
function(alpha){



}
