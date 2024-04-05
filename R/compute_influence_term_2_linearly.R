#' Compute Second Influence Term Assuming Piece-wise Linearity
#'
#' Compute the second term of the influence function. defined as
#' $$
#' \int_{t=a}^b V^{-1}B(t)\hat{E}[Y(t)|\bar{O}(t)_i]dt
#' $$
#'
#' The limits of the integral are inferred from the knots of `base`.
#' The integral is computed on the assumption
#'
#' @param df_i data frame for one individual, should include all observations including baseline.
#' @param expected_value The function to compute the expected value of the outcome model.
#' @param base a [Spline Basis][SplineBasis] object.
#' @inheritDotParams ... passed to impute_patient_df
#'
#' @return
#' A data frame with columns `alpha`, and `influence_term_2`.
#' The latter a list column with each element a vector.
#'
#' @export
#'
#' @examples
compute_influence_term_2_linearly <-
function(
    df_i,
    expected_value,
    base,
    variables,
    ...
){
    assert_that(
        is.data.frame(df_i),
        is(base, 'SplineBasis')
    )


    a <- min(base@knots)
    b <- max(base@knots)
    obase <- orthogonalize(base)

    B1 <- integrate(obase)
    B2 <- integrate(B1)


    patient.df <- df_i |>
        filter(
            !!a <= !!variables$time,
            !!variables$time <= !!b
        )

    times <- c(a, pull(patient.df, variables$time), b)

    left <- impute_patient_df(head(times, -1), df_i, variables = variables, ..., right = FALSE) |> expected_value(alpha)
    right <- impute_patient_df(tail(times, -1), df_i, variables = variables, ...) |> expected_value(alpha)

    dt <- diff(times)

    C1 = (right-left)/dt
    C0 = right-C1*tail(times, -1)
    if(F){

        geom_abline(data=tibble(
            c0 = C0[,2],
            c1 = C1[,2],
            period = factor(seq.int(nrow(C1)))
        ), aes(slope=c1, intercept=c0, col=period, group=period))

        library(ggplot2)


    }

    eB1 <- evaluate(B1, times)
    eB2 <- evaluate(B2, times)

    influence_term_2 <- map(seq_along(alpha), function(i){
        pmap(
            list(
                c0 = C0[,i],
                c1 = C1[,i],
                t1 = head(times, -1),
                t2 = tail(times, -1)
            ),
            function(c0, c1, t1, t2){
                (c0+c1*t2)*evaluate(B1,t2) - (c0+c1*t1)*evaluate(B1,t1) -
                    c1*(evaluate(B2, t2) - evaluate(B2, t1))
            }) |> reduce(`+`)
    })

    tibble(alpha, term2 = influence_term_2)
}
