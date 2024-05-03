estimate_baseline_intensity <-
function(
    intensity.model,
    data = NULL,
    bandwidth = NULL,
    kernel = \(x) 0.75*(1 - (x)**2) * (abs(x) < 1),
    variables
){

    mf <- if(is.null(data)){
        model.frame(intensity.model)
    } else {
        model.frame(terms(intensity.model), data = data)
    }

    new.time <- mf[[1]][,2]
    new.strata <- mf[[attr(terms(intensity.model), 'specials')$strata]] |> as.character()


    ############  Estimated baseline intensities for each stratum ############
    tmp <- terms(intensity.model)
    as.list(tmp)[-c(attr(tmp, 'response'), attr(tmp, 'specials')$strata)]


    data_surv <- survfit(intensity.model, newdata = tibble(!!variables$prev_outcome:=0))
    strata <- data_surv$strata

    # match(mf[[2]], names(strata), nomatch = 1L)

    # Andrew [@halpo]: change this part for a general version (use that general v)

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

    colSums(kernel(time.diff/bandwidth) * strata.equal * cumhaz.data$hazard)/bandwidth
}

if(F){
    estimate_baseline_intensity(object$intensity_model, df_i, kernel = dnorm)


}

