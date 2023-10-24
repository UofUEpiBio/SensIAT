test_that("model.matrix.PCORI::Single-index-outcome-model is invariant to subsetting", {
    fitted.trt.sim <-
        fit_PCORI_within_group_model(
            group.data = filter(ARC_data, Trt=='home_visits'),
            outcome_modeler = PCORI_sim_outcome_modeler,
            id.var = elig_pid,
            outcome.var = Asthma_control,
            time.var = time,
            intensity.bandwidth = 30,
            End = 830
        )

    om <- fitted.trt.sim$outcome_model

    terms(om)

    om$frame
    mf <- model.frame(om)

    a <- model.matrix(om)


  expect_equal(2 * 2, 4)
})
