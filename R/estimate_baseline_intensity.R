globalVariables('cumhaz')

estimate_baseline_intensity <-
function(
    intensity.model,
    data = NULL,
    bandwidth = attr(intensity.model, 'bandwidth'),
    kernel = attr(intensity.model, 'kernel') %||% \(x) 0.75*(1 - (x)**2) * (abs(x) < 1)
){

    mf <- if(is.null(data)){
        model.frame(intensity.model)
    } else {
        model.frame(terms(intensity.model), data = data)
    }

    new.time <- mf[[1]][,2]
    new.strata <- mf[[attr(terms(intensity.model), 'specials')$strata]] |> as.character()


    ############  Estimated baseline intensities for each stratum ############
    data_surv <- survfit(intensity.model, newdata = construct_baseline_df(intensity.model))
    strata <- data_surv$strata

    # the following is to find the baseline intensity function for each strata k
    cumhaz.data <-
        tibble(time = data_surv$time,
               cumhaz = data_surv$cumhaz,
               strata = rep(names(strata), strata)
        ) |>
        group_by(strata) |>
        mutate(hazard = cumhaz - dplyr::lag(cumhaz, default=0, order_by = time))

    if(is.null(bandwidth))
        bandwidth <- dpill(cumhaz.data$time, cumhaz.data$hazard)

    time.diff <- outer(cumhaz.data$time, new.time, `-`)
    strata.equal <- outer(cumhaz.data$strata, new.strata, `==`)

    list(
        baseline_intensity = colSums(kernel(time.diff/bandwidth) * strata.equal * cumhaz.data$hazard)/bandwidth,
        bandwidth = bandwidth,
        kernel = kernel
    )
}

construct_baseline_df <- function(intensity.model){
    df <- model.frame(intensity.model)[1,]
    mm <- model.matrix(intensity.model, data = df)
    mutate(df, across(any_of(colnames(mm)), ~0))
}
