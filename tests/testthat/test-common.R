test_that("pcoriaccel_inner_prod works", {
    pcoriaccel_inner_prod(1:3, 3:1) |> expect_equal(10)
})
test_that("pcoriaccel_outer_sum works", {
    pcoriaccel_outer_sum( c(1,2,3,4,5), c(2,4,6) ) |>
        expect_equal(matrix(c(3,5,7,
                              4,6,8,
                              5,7,9,
                              6,8,10,
                              7,9,11), nrow=5, byrow=TRUE))
})
test_that("pcoriaccel_outer_prod works", {
    pcoriaccel_outer_prod( c(1,2,3,4,5), c(2,4,6) ) |>
        expect_equal(matrix(c(2,4,6,
                              4,8,12,
                              6,12,18,
                              8,16,24,
                              10,20,30), nrow=5, byrow=TRUE))
})
