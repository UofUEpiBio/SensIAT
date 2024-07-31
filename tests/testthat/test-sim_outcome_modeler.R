test_that("model.matrix.PCORI::Single-index-outcome-model is invariant to subsetting", {
    object <-
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
    expect_s3_class(object, 'PCORI_within_group_model')
    expect_contains(
        names(object),
        c('intensity.model', 'outcome.model',
          'data', 'variables', 'End', 'influence', 'coefficients',
          'coefficient.variance', 'base'))
    expect_s3_class(object$outcome.model, 'PCORI::Single-index-outcome-model')
    expect_s3_class(object$intensity.model, 'coxph')
    expect_true(is.list(object$influence))
    expect_true(is.list(object$influence[[1]]))
    expect_contains(names(object$influence[[1]]), c('id', 'alpha', 'term1', 'term2'))
})
