estimate_spline_covariance <-
function(
    spline_fn,
    interval,
    resolution = diff(interval),
    ...
){

    x <- seq(min(interval), max(interval), length.out=resolution)
    d <- abs(diff(interval))/resolution

    B <- t(sapply(x, spline_fn, ...))

    d*(crossprod(head(B, -1)) + crossprod(tail(B,-1)))/2
}
if(F){
    a <- estimate_spline_covariance(spline_fn, c(60, 460), knots = c(59,59,59,59,260,461,461,461,461))
    b <- estimate_spline_covariance(spline_fn, c(60, 460), 1000, knots = c(59,59,59,59,260,461,461,461,461))
    c <- estimate_spline_covariance(spline_fn, c(60, 460), 1e5, knots = c(59,59,59,59,260,461,461,461,461))

    a-b
    b-c
    all.equal(a,b)
    all.equal(b,c)


}

