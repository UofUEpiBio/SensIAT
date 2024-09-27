cumuSIR_new <- function(X, Y, eps = 1e-7)
{
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    number_n <- dim(X)[1]
    number_p <- dim(X)[2]

    Y.CP <- matrix(
        Y[rep(1:number_n, times = number_n), ]<=
            Y[rep(1:number_n, each = number_n), ],
        nrow = number_n, ncol = number_n
    )

    # centralizing covariates
    X.cs <- t(t(X)-colMeans(X))

    # calculating m(y)=\E[X_i 1(Y_i\leq y)]
    m.y <- t(X.cs) %*% Y.CP/number_n
    # calculating K=\E[m(Y_i)m(Y_i)^T]
    Km <- m.y %*% t(m.y)/number_n

    Bhat <- eigen(solve(var(X) + eps*diag(number_p), Km))$vectors

    return(list(basis = Bhat))
}
