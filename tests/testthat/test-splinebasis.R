
test_that("Compiled version of evalute works as expected.", {

    basis <- SplineBasis(c(60,60,60,60,260,460,460,460,460))
    vals <- c( 59.0,60.0,61.0, 259.0,260.0,261.0, 459.0,460.0,461.0 )

    expect_equal(pcoriaccel_evaluate_basis(basis, 59.0), rep(NA_real_, 5))
    expect_equal(pcoriaccel_evaluate_basis(basis, 60.0), c(1, 0, 0, 0, 0))
    expect_equal(pcoriaccel_evaluate_basis(basis, 61.0), evaluate(basis, 61.0)[1,])
    expect_equal(pcoriaccel_evaluate_basis(basis, 259.0), evaluate(basis, 259.0)[1,])
    expect_equal(pcoriaccel_evaluate_basis(basis, 260.0), evaluate(basis, 260.0)[1,])
    expect_equal(pcoriaccel_evaluate_basis(basis, 261.0), evaluate(basis, 261.0)[1,])
    expect_equal(pcoriaccel_evaluate_basis(basis, 459.0), evaluate(basis, 459.0)[1,])
    expect_equal(pcoriaccel_evaluate_basis(basis, 460.0), c(0, 0, 0, 0, 1))
    expect_equal(pcoriaccel_evaluate_basis(basis, 461.0), rep(NA_real_, 5))

})
