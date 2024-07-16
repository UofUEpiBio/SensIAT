test_that("compute_influence_term_2_quadv_cpp_sim is faster but with equal results", {
    time.quadv <- system.time(
        fitted.trt.sim.quadv <-
            fit_PCORI_within_group_model(
                group.data = PCORI_example_data,
                outcome_modeler = PCORI_sim_outcome_modeler,
                id.var = Subject_ID,
                outcome.var = Outcome,
                time.var = Time,
                control = pcori_control(
                    intensity.bandwidth = 30,
                    integration.method = 'quadv',
                    tol = .Machine$double.eps^(1/2)
                ),
                intensity.bandwidth = 30,
                knots = c(60,60,60,60,260,460,460,460,460),
                End = 830
            )
    )

    object <- fitted.trt.sim.quadv

    value <- compute_influence_term_2_quadv_cpp_sim(
        object$data |> filter(Subject_ID == 1),
        alpha = 0.5,
        object$outcome.model,
        object$base,
        object$variables,
        object$outcome.model.centering
    )
    value_prev <- compute_influence_term_2_quadv_sim(
        object$data |> filter(Subject_ID == 1),
        alpha = 0.5,
        object$outcome.model,
        object$base,
        object$variables,
        object$outcome.model.centering
    )

    expect_equal(value, value_prev, ignore_attr = TRUE)

    time.cpp <- system.time(
        fitted.trt.sim.quadvcpp <-
            fit_PCORI_within_group_model(
                group.data = PCORI_example_data,
                outcome_modeler = PCORI_sim_outcome_modeler,
                id.var = Subject_ID,
                outcome.var = Outcome,
                time.var = Time,
                control = pcori_control(
                    intensity.bandwidth = 30,
                    integration.method = 'quadvcpp',
                    tol = .Machine$double.eps^(1/2)
                ),
                intensity.bandwidth = 30,
                knots = c(60,60,60,60,260,460,460,460,460),
                End = 830
            )
    )

    expect_true(all.equal(
        fitted.trt.sim.quadvcpp$influence |> pull(term2) |> do.call(rbind, args=_),
        fitted.trt.sim.quadv   $influence |> pull(term2) |> do.call(rbind, args=_),
        tolerance = .Machine$double.eps^(1/4)
    ))

    skip_on_cran()
    skip_on_ci()
    expect_true(time.cpp[1] < time.quadv[1])
})
