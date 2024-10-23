test_that("Compute Influence term2 fixed point integration", {
    fixed.object <-
        fit_PCORI_within_group_model(
            group.data = PCORI_example_data,
            outcome_modeler = PCORI_sim_outcome_modeler,
            alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
            id.var = Subject_ID,
            outcome.var = Outcome,
            time.var = Time,
            intensity.bandwidth = 30,
            knots = c(60,60,60,60,260,460,460,460,460),
            End = 830,
            influence.args = list(method = 'fixed', delta = 1)
        )
    expect_identical(attr(fixed.object$influence[[3]]$term2, 'fcnt'),
              vector(mode='list', length=100))
    expect_identical(attr(fixed.object$influence[[3]]$term2, 'estim.prec'),
              vector(mode='list', length=100))
})
