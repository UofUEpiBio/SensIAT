
#' Format the original ARC data
#'
#' @param df cleaned ARC data
#' @param id_var A string indicating the patient/subject identifier variable, must be present in `df`.
#' @param treatment_var  A string indicating the variable containing the treatment group, , must be present in `df`.
#' @param visitnumber_var A string indicating an numeric(integer) variable from `df` to identify the visit number within patient.
#' @param outcome_var A string indicating the outcome variable from `df`.
#' @param time_var A string indicating the time variable from `df`.
#' @param spline_seq passed to [spline_fn] to create the spline coefficient matrix.
#' @param V TODO[@yujinggao24]: please clarify.
#' @param knots passed on in output.
#' @param p The number of rows of data.  TODO[@yujinggao24]: Please Verify.
#' @param last_day Last time point of observation, defaults to the maximum time present in `df`.
#'
#' @description
#' This performs additional formatting on the cleaned ARC data.
#'
#' @return
#' A list of three [data frames](data.frame) with class `arc`:
#'
#' * baseline visits,
#' * post-baseline visits, and
#' * a [data frame](data.frame) for use with the [`survival`](package:survival)
#'   package, which has rows indicationg censoring for subjects with fewer than
#'   4 post-baseline assessments.
#'
#' Also attached are attributes:
#'
#'  * spline_seq, the same as input;
#'  * End, equal to input `last_day`;
#'  * V_inverse, $V^{-1}$ where $V=B(t)B(t)'$ and $B(t)$ is the spline matrix created by
#'    applying `spline_fn` to `spline_seq`;
#'  * Weights_term2, $V^{-1}B(t)$;
#'  * knots, passed through from input.
#'  * p, the number of rows of $B(t)$.
#'
#'
#' @export
#'
#' @examples
#' data(ARC_data, package='pcoriRPackage')
#' ARC_formatted <- formatting_fn(
#'     df              = ARC_data,
#'     id_var          = "elig_pid",
#'     treatment_var   = "Trt",
#'     outcome_var     = "Asthma_control",
#'     visitnumber_var = "Visit_number",
#'     time_var        = "time",
#'     V               = 5,
#'     knots           = c(59,59,59,59,260,461,461,461,461),
#'     spline_seq      = 60:460,
#'     last_day = 830
#' )

formatting_fn <- function(
        df,
        id_var,
        treatment_var,
        visitnumber_var,
        outcome_var,
        time_var,
        spline_seq,
        V,
        knots,
        p = length(knots)-4,
        last_day = max(df[[time_var]], na.rm = TRUE)
){

    # Subsetting to the last day and learning the dims
    df <- subset(df, df[[time_var]] <= last_day)
    K  <- dim(df)[1]

    # Lagging
    Prev_outcome <- c(NA_real_, df[[outcome_var]][-1])
    Prev_time    <- c(NA_real_, df[[time_var]][-1])

    # Correcting for not matching
    Prev_outcome[ df[[visitnumber_var]] == 0] <- NA_real_
    Prev_time[ df[[visitnumber_var]] == 0]    <- NA_real_

    df$Prev_outcome <- Prev_outcome
    df$Prev_time    <- Prev_time
    df$Lag_time     <- df[[ time_var ]] - Prev_time
    df$Event        <- 1

    ####  Each subject's baseline visit
    Baseline_df <- subset(df, df[[ visitnumber_var ]] == 0) |>
        subset(select = c(-Event,-Prev_outcome,-Prev_time,-Lag_time))

    ####  Each subject's post-baseline visits
    Visits_df <- subset(df, df[[ visitnumber_var ]] != 0) |>
        subset(select = c(-Event))

    Survival_df <- df

    u <- sort(unique(Survival_df[[ id_var ]]))
    H <- length(u)

    # [2022-12-19] George Seems to be completing the number of observations.
    less_than_V <- as.integer(names(which(
        table(ARC_data[[ id_var ]]) < V
    )))

    for(h in less_than_V){

        # Counting the number of rows
        df_h <- subset(Survival_df, Survival_df[[ id_var ]] == h)
        v    <- dim(df_h)[1] # Number of rows in the obs

        # If less than 5 observations, then repeat the last one
        # for some reason I don't understand
        temp_df <- data.frame(
            elig_pid       = h,
            Asthma_control = NA,
            time           = last_day,
            Trt            = df_h[[treatment_var]][v],
            Visit_number   = v,
            Prev_outcome   = df_h[[ outcome_var ]][v],
            Prev_time      = df_h[[ time_var ]][v],
            Lag_time       = last_day - df_h[[ time_var ]][v],
            Event          = 0
        )

        Survival_df <- rbind(Survival_df, temp_df)

    }

    Survival_df <- Survival_df[order(Survival_df[[ id_var ]]),]
    Survival_df <- subset(Survival_df, Survival_df[[ visitnumber_var ]] != 0)

    ##### the p x 401 matrix with columns B(t), t=60,...,460
    B_t_matrix <- matrix(
        sapply(
            spline_seq,
            spline_fn,
            knots = knots
        ),
        byrow = FALSE,
        nrow  = p
    )

    #### we approximated the integral V=\int_t B(t)B(t)'dt using sums of rectangles of width 1
    V         <- crossprod(t(B_t_matrix))
    V_inverse <- solve(V)

    Weights_term2 <- V_inverse %*% B_t_matrix

    structure(
        list(BaseLine = Baseline_df,
             Visits   = Visits_df,
             Survival = Survival_df
        ),
        spline_seq    = spline_seq,
        End           = last_day,
        V_inverse     = V_inverse,
        Weights_term2 = Weights_term2,
        knots         = knots,
        p             = p,
        class         = "arc"
    )

}
