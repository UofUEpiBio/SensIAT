test_that("outcome model structure", {
    runif(1)
    old.seed <- .Random.seed
    source(system.file("examples/basic.R", package='SensIAT'))
    new.seed <- .Random.seed
    expect_identical(old.seed, new.seed)

    expect_s3_class(object, 'SensIAT_within_group_model')
    expect_contains(
        names(object),
        c('models', 'data', 'variables', 'End', 'influence',
          'alpha', 'coefficients', 'coefficient.variance',
          'args', 'base', 'V_inverse'))
    expect_true(is.list(object$models))
    expect_s3_class(object$data, 'data.frame')
    expect_true(is.list(object$variables) && rlang::is_named(object$variables))

    expect_true(is.list(object$influence))
    expect_length(object$influence, length(object$alpha))
    expect_true(is.list(object$influence[[1]]))
    expect_named(object$influence[[1]], c('term1', 'term2', 'id', 'alpha'), ignore.order=TRUE)

    expect_true(is.list(object$models))
    expect_named(object$models, c('intensity', 'outcome'), ignore.order=TRUE)
    expect_s3_class(object$models$outcome, 'SensIAT::Single-index-outcome-model')
    expect_s3_class(object$models$intensity, 'coxph')

    expect_equal(object$alpha, 0)
    expect_s4_class(object$base, 'SplineBasis')

    expect_true(is.list(object$args))
    expect_named(object$args, c('intensity', 'outcome', 'influence'), ignore.order=TRUE)

    expect_true(is.list(object$args$intensity))
    expect_true(is.list(object$args$outcome))
    expect_true(is.list(object$args$influence))

})
