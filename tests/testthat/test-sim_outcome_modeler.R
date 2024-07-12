test_that("model.matrix.PCORI::Single-index-outcome-model is invariant to subsetting", {
    fitted.trt.sim <-
        fit_PCORI_within_group_model(
            group.data = PCORI_example_data,
            outcome_modeler = PCORI_sim_outcome_modeler,
            id.var = Subject_ID,
            outcome.var = Outcome,
            time.var = Time,
            intensity.bandwidth = 30,
            knots = c(60,60,60,60,260,460,460,460,460),
            End = 830
        )
    expect_s3_class(fitted.trt.sim, 'PCORI_within_group_model')
    expect_contains(
        names(fitted.trt.sim),
        c('intensity.model', 'outcome.model', 'outcome.model.centering',
          'data', 'variables', 'End', 'influence', 'coefficients',
          'coefficient.variance', 'control', 'base'))
    expect_is(fitted.trt.sim$outcome_model, 'PCORI::Single-index-outcome-model')

})
