SensIAT_data_transform <-
    function(data,
             id, outcome, time,
             end.time = pull(summarize(data, max({{time}}, na.rm = TRUE)))
             ){
        id.var <- ensym(id)
        outcome.var <- ensym(outcome)
        time.var <- ensym(time)
        vars <- list(
            outcome = outcome.var,
            id = id.var,
            time = time.var
        )
        force(end.time)
        if(is.null(end.time)) {
            end.time <- pull(summarize(data, max({{time}}, na.rm = TRUE)))
        }
        data |>
            filter({{time}} <= end.time) |>
            rename(
                ..id.. = !!id.var,
                ..time.. = !!time.var,
                ..outcome.. = !!outcome.var
            ) |>
            group_by(..id..) |>
            mutate(
                ..visit_number.. = seq_along(..time..)
            ) |>
            ungroup() |>
            complete(..id.., ..visit_number..,
                     fill = list(..time.. = end.time,..outcome.. = NA_real_)
            ) |>
            group_by(..id..) |>
            arrange(..id.., ..visit_number..) |>
            mutate(
                ..time..            := as.double(..time..),
                ..prev_outcome..    := lag(..outcome.., order_by = ..time..),
                ..prev_time..       := lag(..time.., order_by =  ..time.., default = 0),
                ..delta_time..      := ..time.. - lag(..time.., order_by =  ..time.., default = 0)
            ) |>
            ungroup() |>
            structure(variables = vars)
    }
