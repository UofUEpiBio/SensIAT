compute_slope <-
    function(outcome.model,
             time.vars,
             start = NULL,
             delta = 1,
             ... #< IGNORED
    ){
        vars <- all.vars(delete.response(terms(outcome.model)))
        if(is.null(start)) {
            # Default start values
            start <- rep(list(0), length(vars))
            names(start) <- vars
            start <- as.data.frame(start)
        } else {
            if (is.numeric(start)){
                if (length(start) == 1){
                    start <- rep(start, length(vars))
                    names(start) <- vars
                } else {
                    if(rlang::is_named2(start))
                        rlang::abort("start vector of length > 1 must be named.")
                    start <- as.list(start)
                }
            } else
                if (is.list(start)) {
                    if (rlang::is_named2(start)) {
                        start <- as.data.frame(start)[names(start) %in% vars]
                    } else {
                        rlang::abort("start list must be named.")
                    }
                    start <- as.data.frame(start)[names(start) %in% vars]
                } else
                    if(is.data.frame(start)) {
                        start <- start[vars]
                    } else {
                        stop("start must be a numeric vector, list, data frame, or NULL")
                    }
        }
        end = start
        for(var in names(start)){
            if (var %in% time.vars)
                end[[var]] <- start[[var]] + delta
            else if (is.na(start[[var]]))
                end[[var]] <- start[[var]] <- 0
        }
        a <- model.matrix(delete.response(terms(outcome.model)), data = start)
        b <- model.matrix(delete.response(terms(outcome.model)), data = end)

        as.vector((b-a)/delta)
    }
