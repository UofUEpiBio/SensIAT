#' Produce fitted model for group (treatment or control)
#'
#' Produces a fitted model that may be used to produce estimates of mean and
#' variance for the given group.
#'
#' This function should be agnostic to whether it is being provided a
#' treatment or control group.
#'
#' @param group.data The data for the group that is being analyzed.
#'          Preferably passed in as a single tibble that internally is
#'          subsetted/filtered as needed.
#' @param End The end time for this data analysis, we need to set the default value as the
#'           max value of the time
#' @param kernel
#' @param nmk
#' @param outcome_modeler A separate function that may be swapped out to switch
#'          between negative-binomial, single index model, or another we will
#'          dream up in the future.
#' @param ... add parameters as needed or use this to pass forward into the
#'          outcome_modeler.
#'
#' @return
#'      Should return everything needed to define the fit of the model.
#'      This can then be used for producing the estimates of mean, variance,
#'      and in turn treatment effect.
#'
#' @export
#'
#' @examples
#'
#' fitted.trt.nb <-
#'     fit_PCORI_within_group_model(
#'         group.data = filter(ARC_data, Trt=='home_visits'),
#'         id.var = elig_pid,
#'         outcome.var = Asthma_control,
#'         time.var = time,
#'         intensity.bandwidth = 30,
#'         End = 830,
#'         outcome_modeler = glm.nb
#'     )
#'
#' fitted.trt.sim <-
#'     fit_PCORI_within_group_model(
#'         group.data = filter(ARC_data, Trt=='home_visits'),
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         id.var = elig_pid,
#'         outcome.var = Asthma_control,
#'         time.var = time,
#'         intensity.bandwidth = 30,
#'         End = 830
#'     )
#'
fit_PCORI_within_group_model <- function(
        group.data,
        outcome_modeler,
        id.var,
        outcome.var,
        time.var,
        intensity.covariates = ~1,
        outcome.covariates = ~-1,
        End = max({{time.var}}, na.rm = TRUE) + 1,
        ...
){
    id.var <- ensym(id.var)
    outcome.var <- ensym(outcome.var)
    time.var <- ensym(time.var)

    outcome_modeler <- match.fun(outcome_modeler)
    End <- rlang::enexpr(End)
    End <- rlang::eval_tidy(End, data = group.data, env =parent.frame())

    vars <- list(
        id = id.var,
        time = time.var,
        outcome = outcome.var,
        prev_outcome = rlang::sym(glue::glue("lag({outcome.var})")),
        prev_time = rlang::sym(glue::glue("lag({time.var})")),
        delta_time = rlang::sym(glue::glue("Δ({time.var})")),
        norm_time = rlang::sym(glue::glue("scale({time.var})")),
        norm_delta_time = rlang::sym(glue::glue("scale(Δ({time.var}))"))
    )

    # For example, if we want to fit models for the treatment group
    # Then the group.data should includes parameter "End"
    # Dataframe: "data_baseline_hv", "data_visits_hv", "data_survival_hv"

    group.data2 <- filter(group.data, !!time.var <= End)

    # data_formatted <- formatting_fn(df = group.data2,
    #      id_var          = "elig_pid",
    #      treatment_var   = "Trt",
    #      outcome_var     = "Asthma_control",
    #      visitnumber_var = "Visit_number",
    #      time_var        = "time",
    #      V               = 5,
    #      knots           = c(59,59,59,59,260,461,461,461,461),
    #      spline_seq      = 60:460,
    #      last_day = 830
    #  )
    # data_baseline_hv <- data_formatted[[1]]
    # data_visits_hv <- data_formatted[[2]]
    # data_survival_hv <- data_formatted[[3]]

    u_hv <- group.data2 |> select(!!id.var) |> distinct() |> pull()
    N <- pull(summarize(group.data2, n_distinct(!!id.var)))
    # N <- dim(data_baseline_hv)[1]


    # Andrew [@halpo]: the user can only change this part, "Prev_outcome+strata(Visit_number)"

    ######   Andersen-Gill model stratifying by assessment number

    model.data <-
        rlang::inject(!!outcome.var ~ !!id.var + !!time.var + !!rlang::f_rhs(intensity.covariates)) |>
        model.frame(data=filter(group.data, (!!time.var) < !!End)) |>
        arrange(!!id.var, !!time.var) |>
        group_by(!!id.var) |>
        mutate(
            visit.number = seq_along(!!time.var)
        ) |>
        ungroup() |>
        complete(!!id.var, visit.number, fill = tibble::lst(!!time.var := !!End)) |>
        group_by(!!id.var) |>
        arrange(!!id.var, visit.number) |>
        mutate(
            !!vars$prev_outcome    := lag(!!outcome.var, order_by = !!time.var),
            !!vars$prev_time       := lag(!!time.var, order_by =  !!time.var, default = 0),
            !!vars$delta_time      := !!time.var - lag(!!time.var, order_by =  !!time.var, default = 0),
            across(all_of(time.var), coalesce, !!End)
        ) |>
        ungroup() |>
        filter(!(visit.number > 1 & is.na(lag(!!outcome.var))))

    centering.statistics <-
        summarize( ungroup(filter(model.data, !!time.var > 0, !is.na(!!outcome.var))),
            "mean({time.var})" := mean(!!time.var),
            "sd({time.var})" := sd(!!time.var),
            "mean({vars$delta_time}))" := mean(!!vars$delta_time),
            "sd({vars$delta_time}))" := sd(!!vars$delta_time)
        )

    model.data <- model.data |>
        mutate(
            !!vars$norm_time := (!!vars$time - !!pull(centering.statistics, 1))/!!pull(centering.statistics, 2),
            !!vars$norm_delta_time := (!!vars$delta_time - !!pull(centering.statistics, 3))/!!pull(centering.statistics, 4)
        )
    ###########################   Model 1: Intensity model  ######################
    intensity.model <-
        rlang::inject(coxph(
            Surv(!!vars$prev_time,
                 !!time.var,
                 !is.na(!!outcome.var))~!!vars$prev_outcome+strata(visit.number),
            id = !!id.var,
            data = filter(model.data, !!time.var > 0)
        ))

    # gamma <- Int_model$coefficients # parameter lambda in lambda(t, O(t))
    # we need this gamma value


    ############  Estimated baseline intensities for each stratum ############
    if(F){
        # data_surv <- survfit(intensity.model, newdata=data.frame(Prev_outcome=0))
        data_surv <- survfit(intensity.model, newdata=tibble(!!vars$prev_outcome := 0))
        strata <- data_surv$strata

        # Andrew [@halpo]: change this part for a general version (use that general v)

        # the following is to find the baseline intensity function for each strata k

        cumhaz.data <-
            tibble(time = data_surv$time,
                   cumhaz = data_surv$cumhaz,
                   strata = factor(rep(names(strata), strata), levels = names(strata))
                   ) |>
            group_by(strata) |>
            mutate(hazard = cumhaz - lag(cumhaz, default=0, order_by = time))
        base_intens1 <- cumhaz.data |>
            group_by(strata) |>
            group_map(rlang::inject(~purrr::map_dbl(seq.int(!!End), lambda0_fn, b=intensity.bandwidth, surv = .x)))

        tmp <- outer(data_surv$time, seq.int(End), `-`)
        increment <- cumhaz.data$hazard

        base_intens <- apply(
            0.75*(1 - (tmp/intensity.bandwidth)**2) *
                (abs(tmp) < intensity.bandwidth) *
                cumhaz.data$hazard,
            2, tapply, cumhaz.data$strata, sum
        )

        # v1=strata[1]
        # cumhaz_v1=data.frame(time <- data_surv$time[1:v1],
        #                      cumhaz <- data_surv$cumhaz[1:v1])
        # base_intens_v1=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v1)
        #
        # v2=strata[2]
        # cumhaz_v2=data.frame(time=data_surv$time[(v1+1):(v1+v2)],
        #                      cumhaz=data_surv$cumhaz[(v1+1):(v1+v2)])
        # base_intens_v2=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v2)
        #
        # v3=strata[3]
        # cumhaz_v3=data.frame(time=data_surv$time[(v1+v2+1):(v1+v2+v3)],
        #                      cumhaz=data_surv$cumhaz[(v1+v2+1):(v1+v2+v3)])
        # base_intens_v3=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v3)
        #
        # v4=strata[4]
        # cumhaz_v4=data.frame(time=data_surv$time[(v1+v2+v3+1):(v1+v2+v3+v4)],
        #                      cumhaz=data_surv$cumhaz[(v1+v2+v3+1):(v1+v2+v3+v4)])
        # base_intens_v4=sapply(1:End,lambda0_fn,b=30,surv=cumhaz_v4)
    }


    #############  Model 2: Outcome model - single index model   ############
    {

        # Andrew [@halpo]: keep the single index model as a flexible model

        # Question: Does this make sense to scale the Lag_time?  should we be lagging scaled time instead?
        outcome.formula <- rlang::inject(
            !!outcome.var~
                ns(!!vars$prev_outcome, df=3) +
                !!vars$norm_time +
                !!vars$norm_delta_time +
                !!rlang::f_rhs(outcome.covariates)
        )
        # data_visits_hv1 <- cbind(data_visits_hv,
        #                          ns1 = ns(data_visits_hv$Prev_outcome, df = 3)[, 1],
        #                          ns2 = ns(data_visits_hv$Prev_outcome, df = 3)[, 2],
        #                          ns3 = ns(data_visits_hv$Prev_outcome, df = 3)[, 3],
        #                          time_scale = scale(data_visits_hv$time),
        #                          Lag_time_scale = scale(data_visits_hv$Lag_time))

        # time_mean <- mean(data_visits_hv1$time)
        # time_sd <- sd(data_visits_hv1$time)
        # Lag_time_mean <- mean(data_visits_hv1$Lag_time)
        # Lag_time_sd <- sd(data_visits_hv1$Lag_time)



        # Xi <- data.frame(data_visits_hv1$ns1,
        #                  data_visits_hv1$ns2,
        #                  data_visits_hv1$ns3,
        #                  data_visits_hv1$time_scale,
        #                  data_visits_hv1$Lag_time_scale)
        #
        # Xi <- as.matrix(Xi)
        # Yi <- as.matrix(data_visits_hv1$Asthma_control, ncol = 1)

        Xi <- model.matrix(outcome.formula, data = filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))
        Yi <- model.response(model.frame(outcome.formula, data = filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))) |>
                as.matrix(ncol=1)

        # Andrew [@halpo]: add these functions into the R package
        # new method to get the SIM's estimation
        # SDR1 <- cumuSIR_new(X = Xi, Y = Yi)
        # Outcome_model <- SIDRnew(X = Xi, Y = Yi, initial = SDR1$basis[, 1], kernel = "dnorm", method = "nmk")
    }
    outcome.model <- outcome_modeler(outcome.formula, data = filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))

    structure(list(
        intensity_model = intensity.model,
        # baseline_intensity = base_intens,
        outcome_model = outcome.model,
        outcome_model_centering = centering.statistics,
        data = model.data,
        variables = vars,
        End = End
    ), class = "PCORI_within_group_model")
}



#' Prediction method for `PCORI_within_group_model` objects
#'
#' @param object a `PCORI_within_group_model` object.
#' @param time The time(s) to evaluate at.
#' @param alpha The values for the sensitivity parameter to evaluate at.
#' @param spline_fn Function to to generate the values for the spline functions
#'      with knots given by `knots`.  Must accept time as the first argument and
#'      knots as the second.
#' @param spline_seq The time points used to approximate the integral in term 2
#'          of the influence function.
#' @param knots The knots used for the spline function.
#'          Includes both boundary and internal knots.
#' @param integration.method Method for integration when computing the second influence term.  See Details.
#' @param integral.resolution the number of points to use for numerical integration.
#'
#'  Evaluate the fitted model, `object`, at each combination of `time` and
#'  `alpha`.
#'
#'
#' @return A tibble/data.frame with the components: `time`, `alpha`, `mean`, `var`.
#' where `time` and `alpha` are the combinations of the respective input(s) and
#' `mean` and `var` are the estimated mean and variance of the response for the
#' given model.
#'
#'  @details
#'  For `integration.method` when computing the integral in term two of the
#'  influence function `linear` approximates the expected value as piece-wise
#'  linear and is much faster, `numerical` uses traditional numerical
#'  integration.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' fitted.trt.nb <-
#'     fit_PCORI_within_group_model(
#'         group.data = filter(ARC_data, Trt=='home_visits'),
#'         id.var = elig_pid,
#'         outcome.var = Asthma_control,
#'         time.var = time,
#'         intensity.bandwidth = 30,
#'         End = 830,
#'         outcome_modeler = glm.nb
#'     )
#' predict(fitted.trt.nb, time = c(90, 180), alpha = c(-0.6, -0.3, 0, 0.3, 0.6)
#'        , spline_fn, spline_seq=seq(60, 460, by=1), knots=c(59,59,59,59,260,461,461,461,461)
#' )
#'
#' fitted.trt.sim <-
#'     fit_PCORI_within_group_model(
#'         group.data = filter(ARC_data, Trt=='home_visits'),
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         id.var = elig_pid,
#'         outcome.var = Asthma_control,
#'         time.var = time,
#'         End = 830
#'     )
#' time.pw <- system.time({
#' pred.pw <- predict(fitted.trt.sim, time = c(90, 180),
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     intensity.bandwidth = 30,
#'     knots=c(60,60,60,60,260,460,460,460,460),
#'     integration.method = 'piecewise'
#' )
#' })
#' time.num <- system.time({
#' pred.num <- predict(fitted.trt.sim, time = c(90, 180),
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     intensity.bandwidth = 30,
#'     knots=c(60,60,60,60,260,460,460,460,460),
#'     integration.method = 'numerical'
#' )
#' })
#' time.lin <- system.time({
#' pred.lin <- predict(
#'     fitted.trt.sim, time = c(90,180),
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     knots=c(60,60,60,60,260,460,460,460,460),
#'     integration.method = 'linear'
#' )
#' })
#' pred.num$Beta_hat[[2]]
#' pred.lin$Beta_hat[[3]]
#' pred.num$term1[[2]] |> head()
#' pred.lin$term1[[3]] |> head()
#' pred.num$term2[[2]] |> head()
#' pred.lin$term2[[3]] |> head()
#'
#'
`predict.PCORI_within_group_model` <-
function(object, time, alpha,
         knots,
         intensity.bandwidth = NULL,
         integral.resolution = 1000,
         integration.method = c("linear", "numerical", "piecewise", 'quadv'),
         ...
         ){
    assert_that(
        is(object, 'PCORI_within_group_model'),
        is.numeric(time),
        is.numeric(alpha),
        is.numeric(knots),
        is.count(integral.resolution)
    )
    integration.method = match.arg(integration.method)
    a <- min(knots)
    b <- max(knots)
    base <- OBasis(knots)
    self_contained_spline <- \(x)evaluate(base, x)

    # Compute value of the influence function: -----------------------------
    influence <- purrr::map_dfr(alpha, function(alpha){
        object$data |>
            group_by(!!object$variables$id) |>
            group_modify(
                compute_influence_for_one_alpha_and_one_patient,
                alpha = alpha,
                object = object,
                base=base,
                integration.method = integration.method,
                numerical.ingtegration.resolution = integral.resolution,
                ...,
                .keep=TRUE
            )
    })


    B_t <- evaluate(base, time)

    # Results
    influence |>
        group_by(alpha) |>
        summarize(
            term1 = list(term1 |> reduce(rbind) |> `rownames<-`(NULL)),
            term2 = list(term2 |> reduce(rbind) |> `rownames<-`(NULL)),
            IF = list(influence |> reduce(rbind) |> `rownames<-`(NULL))
        ) |>
        mutate(
            Beta_hat = purrr::map(IF, colMeans),
            Var_beta = purrr::map2(IF, Beta_hat, ~tcrossprod(t(.x) - .y)/(nrow(.x)^2))
        ) |>
        mutate(
            # Time specific estimates.
            mean_t = purrr::map(Beta_hat, ~B_t %*% .x),
            var_t = purrr::map(Var_beta, ~B_t %*% .x %*% t(B_t))
        )
}
compute_influence_for_one_alpha_and_one_patient <-
function(
    df_i,
    alpha,
    object,
    base,
    integration.method = c('piecewise', 'numerical', 'linear', 'quadv'),
    ...,
    numerical.ingtegration.resolution = 1000
){
    if (getOption('PCORI::do_arg_checks', TRUE))
        assert_that(
            rlang::is_atomic(alpha), is.numeric(alpha),
            is(object, 'PCORI_within_group_model'),
            is.data.frame(df_i),
            is(base, "OrthogonalSplineBasis"),
            assertthat::is.count(numerical.ingtegration.resolution)
        )
    integration.method <- match.arg(integration.method)
    # id <- unique(pull(df_i, !!object$variables$id))
    # if (getOption('PCORI::do_arg_checks', TRUE))
    #     assert_that(rlang::is_atomic(id))
    df_i[!is.na(pull(df_i, object$variables$prev_outcome)), 'baseline_lambda'] <-
        estimate_baseline_intensity(object$intensity_model, df_i[!is.na(pull(df_i, object$variables$prev_outcome)), ])

    df.in.range <- df_i |>
        filter(
            !!min(base@knots) <= !!object$variables$time,
            !!object$variables$time <= !!max(base@knots)
        )
    if(nrow(df.in.range) == 0) return(tibble())
    # baseline_lambda <- suppressWarnings(
    #     estimate_baseline_intensity(object$intensity_model, df_i, intensity.bandwidth)
    # )
    term1 <- df_i |>
        mutate(
            Exp_gamma = exp((!!coef(object$intensity_model))*!!object$variables$prev_outcome),
        ) |>
        filter(
            !!min(base@knots) <= !!object$variables$time,
            !!object$variables$time <= !!max(base@knots)
        ) |>
        pcori_conditional_means(
            object$outcome_model, alpha, new.data = _
        ) |>
        mutate(
            Term1_unweighted =
                (!!(object$variables$outcome)-E_Y_past)/
                (baseline_lambda*Exp_gamma* exp(-alpha*!!(object$variables$outcome))*E_exp_alphaY)
        ) |>
        summarize(
            term1 =
                list(crossprod(evaluate(base, .data$time), .data$Term1_unweighted)),
        ) |>
        pull(term1) |> unlist()
    term2 <-
        if (integration.method == 'piecewise') {
            numerically_integrate_influence_term_2_for_one_alpha_and_one_patient_piecewise(
                df_i, object=object, alpha=alpha, base=base,
                ...
            ) |> as.vector()
        } else if (integration.method == 'numerical') {
            numerically_integrate_influence_term_2_for_one_alpha_and_one_patient(
                df_i, object, alpha, base,
                resolution = numerical.ingtegration.resolution
            ) |> as.vector()
        } else if (integration.method == 'linear') {
            compute_influence_term_2_linearly(df_i, alpha=alpha, object=object, base=base) |>
                pull(term2) |> unlist() |> as.vector()
        } else if (integration.method == 'quadv') {
            compute_influence_term_2_quadv(df_i, alpha=alpha, object=object, base=base)
        }
    influence <- term1 + term2
    tibble(
        alpha,
        term1 = list(term1),
        term2 = list(term2),
        influence = list(influence)
    )
}
if(F){
    compute_influence_for_one_alpha_and_one_patient(
        df_i, 0, object, base = base, integration.method = 'piecewise', resolution.within.period = 250) |>
        pull(term2)
    compute_influence_for_one_alpha_and_one_patient(
        df_i, 0, object, base = base, integration.method = 'numerical') |>
        pull(term2)
    compute_influence_for_one_alpha_and_one_patient(
        df_i, 0, object, base = base, integration.method = 'linear') |>
        pull(term2)
    compute_influence_for_one_alpha_and_one_patient(
        df_i <- object$data |> filter(elig_pid == unique(elig_pid)[6])
        , 0, object, base = base, integration.method = 'linear') |>
        pull(term1)
}


estimate_influence_term_2 <-
function(object, expected_value, alpha, spline, a, b, resolution, var.map, ...){

    eval.times <- seq(a,b,length.out=resolution)
    B <- spline(eval.times)
    V <- crossprod(B)*(b-a)/resolution
    V_inverse <- solve(V)

    time_mean <- object$outcome_model_centering[[1]]
    time_sd   <- object$outcome_model_centering[[2]]
    Δ_time_mean <- object$outcome_model_centering[[3]]
    Δ_time_sd   <- object$outcome_model_centering[[4]]

    for_one <- function(df_i, g_i, ...){
        outcomes <- c( pull(df_i, !!(object$variables$prev_outcome)),
                       pull(df_i, tail(!!(object$variables$outcome),1)))

        spline_df_est <-
            tibble(
                time = eval.times,
                period = as.numeric(cut(eval.times, c(-Inf, pull(df_i, !!object$variables$time), Inf))),
            ) |>
            mutate(
                delta_time := time - c(0, df_i$time)[period],
                norm_time = (time - time_mean)/time_sd,
                norm_delta_time = (delta_time - Δ_time_mean)/Δ_time_sd,
                prev_outcome = (!!outcomes)[period],
                outcome = 0
            ) |>
            rename(any_of(rlang::set_names(names(object$variables), sapply(object$variables, deparse))))

        Ey = expected_value(spline_df_est, alpha = alpha)

        tibble(
            term2 = list(V_inverse %*% (crossprod(B, Ey) -
                crossprod(head(B,1), head(Ey,1))/2 -
                crossprod(tail(B,1), tail(Ey,1))/2) * (b-a)/resolution)
        )
    }

    tmp <- object$data |> group_by(!!var.map$id) |> group_modify(for_one, alpha = alpha)

    return(mutate(tmp, alpha = alpha))
    # return(purrr::reduce(tmp, `+`)*(b-a)/resolution/length(tmp))

    if(F)for(i in 1:N_hv){
        # print(i)
        df_i1 = object$data |>
            filter(visit.number > 1, !!(object$variables$id)==u_hv[i]) |>
            arrange(!!object$variables$time)

        min_time <- min(pull(df_i1, !!object$variables$time))
        max_time <- max(pull(df_i1, !!object$variables$time))



        outcomes <- c( pull(df_i1, !!(object$variables$prev_outcome)),
                       pull(df_i1, tail(!!(object$variables$outcome),1)))


        # Construct estimating data frame for splines.
        # TODO: handling for extra variables in model.

        Time_means_single <-
            map(alpha, pcori_conditional_means,
                       model = object$outcome_model,
                       new.data = spline_df_est
            )

        E_Y_past_mat <- Time_means_single |> map(pull, 'E_Y_past') |>
            do.call(cbind, args=_)

        # Term2_mat[i, ] <- Weights_term2 %*% E_Y_past_mat
        Term2[i,,] <- Weights_term2 %*% E_Y_past_mat

        # Compute Influence Function, subject specific. ----------------------------
        #######   Subject-specific p x 1 influence function  IF(O_i)

        w = which(pull(Visits_df_a, !!object$variables$id) == u_hv[i])
        Visits_df_a[w,'alpha']
        IF[i,,]=colSums(Term1_mat_2[w,,drop=FALSE]) + Term2[i,,]
    }

}


if(F){
    object <- fitted.trt.sim
    alpha <- c(-0.5, 0, 0.5)

    expected_value <- \(data){matrix(
                pull(pcori_conditional_means(object$outcome_model, alpha = alpha, new.data = data), 'E_Y_past'),
                nrow = nrow(data)
            )}




    spline <- \(x)t(sapply(x, spline_fn, knots = c(59,59,59,59,260,461,461,461,461)))

    a <- 60
    b <- 460
    resolution <- 1000

    df_i <- filter(object$data, elig_pid == elig_pid[[1]])
}


