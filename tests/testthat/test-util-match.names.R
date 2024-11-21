test_that("match.names", {
    expect_identical(
        match.names(list(a=1, b=2), c('apple', 'banana')),
        list(apple=1, banana=2)
    )
    expect_identical(
        match.names(list(a=1, b=2, c=3), c('apple', 'banana')),
        list(apple=1, banana=2, c=3)
    )
    expect_identical(
        match.names(list(a=1, b=2, c=3, b=4), c('apple', 'banana')),
        list(apple=1, banana=2, c=3, banana=4)
    )
    expect_warning(
        match.names(list(a=1, b=2, c=3, ba=4), c('apple', 'banana', 'cherry'))
    )
    expect_error(
        match.names(list(a=1, b=2, c=3, ba=4), c('apple', 'banana', 'cherry'), FALSE)
    )
})
