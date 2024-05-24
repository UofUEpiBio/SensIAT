NW_test <- function(Xb, Y, xb, y, h, kernel = "K2_Biweight"){

    # select kernel
    if(kernel == "dnorm"){
        K <- function(x, h){dnorm(x/h, 0, 1)} # Gaussian
    } else if(kernel == "K2_Biweight"){
        K <- function(x, h){15/16*(1-(x/h)^2)^2 * (abs(x) <= h)} # K2_biweight
    } else if(kernel=="K4_Biweight"){
        K <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }# K4_biweight
    }

    # Compute the kernel applied to each pair of Xb_i and x_j
    # equivalent to outer(Xb, xb, \(a,b){K(a-b, h)})
    # Kxb[i,j] = K(Xb[i] - xb[j], h)
    Kxb <- sapply(xb, function(x, Xb) K(Xb-x, h), Xb=Xb)

    # Ylty[i,k] = Y[i] <= y[k]
    Ylty <- sapply(y, function(x, Y) 1*(Y <= x), Y=Y)

    # denom[j] = sum_i K(Xb[i] - xb[j], h)
    denom <- colSums(Kxb)

    # fyxb[j,k] = sum_i K(Xb[i] - xb[j], h) * 1(Y[i] <= y[k]) / sum_i K(Xb[i] - xb[j], h)
    fyxb <- (denom!=0)*crossprod(Kxb, Ylty)/(denom + (denom==0))

    return(fyxb)
}


# Define the test data
set.seed(20240522)
n <- 100
Xb <- rnorm(n)
Y <- round(rnorm(n), 1)
xb <- seq(-2.5, 2.5, 0.5)
y <- sort(unique(Y))

# parameters
n_y <- length(y)
n_xb <- length(xb)

# Compute the estimated empirical distribution function for the given xb
Fhat <- NW_test(Xb, Y, xb, y, h=1)

# sanity checks
# Check correct dimensions
all(dim(Fhat) == c(n_xb, n_y))
# Check for monotonic increasing
all(apply(Fhat, 1, \(f){all(head(f, -1) <= tail(f, -1))}))
# check all values are between 0 and 1
all(Fhat >= 0)
#all(Fhat <= 1) # failure
all(Fhat - 1 <= .Machine$double.eps) #< works


# save data for results
save(Fhat, file = "../Fhat_test.RData")



