

fit <- survival::survfit(Surv(time, status) ~ x, data = survival::aml)

expect_equal( lambda0_fn(15, b=10, fit), 0.0446332792207792)
expect_error( lambda0_fn(15:45, b=10, fit))

# either allow this
expect_error( lambda0_fn(fit$time, b=10, fit))
# Or allow this
expect_equal(length(lambda0_fn(fit$time, b=10, fit)), length(fit$time))
# but not both


expect_identical( lambda0_fn(15, b=0, fit), NaN)
expect_identical( lambda0_fn(15, b=Inf, fit), 0)


# Trivial Constructed example

faux.surv <-
    list(
        time = 0:100,
        cumhaz = seq(0,3, length.out = 101)
    )

expect_equal(lambda0_fn(  0, 5, faux.surv), 0.0297) #< What should this equal at the boundary of the coverage?
expect_equal(lambda0_fn(  5, 5, faux.surv), 0.0297)
expect_equal(lambda0_fn( 25, 5, faux.surv), 0.0297)
expect_equal(lambda0_fn( 50, 5, faux.surv), 0.0297)
expect_equal(lambda0_fn( 75, 5, faux.surv), 0.0297)
expect_equal(lambda0_fn( 95, 5, faux.surv), 0.0297)
expect_equal(lambda0_fn(100, 5, faux.surv), 0.0297) #< What should this equal at the boundary of the coverage?

# Should the smooth estimate vary with bandwidth when hazard has constant slope?
lambda0_fn(  50, 10, faux.surv)



lambda0_fn(300, 5, faux.surv) #< Should this be error, zero, NA, NaN?





