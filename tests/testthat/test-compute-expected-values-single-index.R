# Validation tests for compute_SensIAT_expected_values.SensIAT::Single-index-outcome-model
#
# Strategy: compare the optimized path (cached Xb_train + pcoriaccel_NW_expectations)
# against a reference computed from the validated pcoriaccel_NW CDF function.
# pcoriaccel_NW itself is tested against the pure-R NW_test reference in test-NW.R.

# Pure-R reference for E[exp(a*Y)|xb] and E[Y*exp(a*Y)|xb] derived from NW CDF matrix
ref_expectations <- function(Xb_train, Yi, xb_query, alpha, h, kernel = "K2_Biweight") {
    y_seq  <- sort(unique(Yi))
    Fhat   <- pcoriaccel_NW(Xb_train, Yi, xb_query, y_seq, h, kernel)  # Q x M
    # PMF for each query point j: P(Y = y_seq[m] | xb[j])
    # Fhat[j, m] = P(Y <= y_seq[m] | xb[j]), so pmf[j, m] = Fhat[j,m] - Fhat[j,m-1]
    pmf <- cbind(Fhat[, 1L, drop = FALSE],
                 Fhat[, -1L, drop = FALSE] - Fhat[, -ncol(Fhat), drop = FALSE])
    exp_a_y   <- exp(alpha * y_seq)          # M-vector
    y_exp_a_y <- y_seq * exp_a_y             # M-vector
    E_exp  <- as.vector(pmf %*% exp_a_y)     # Q-vector
    E_Yexp <- as.vector(pmf %*% y_exp_a_y)   # Q-vector
    list(E_exp = E_exp, E_Yexp = E_Yexp)
}

# Fit a small single-index outcome model for testing.
# initial is supplied explicitly to avoid MAVE::mave.compute being called on
# a tiny dataset (it can return NA, causing `if (initial[1] < 0)` to error).
make_test_model <- function(seed = 1L, n = 60L) {
    set.seed(seed)
    t_vec <- rep(seq(0, 1, length.out = 4L), n)
    subj  <- rep(seq_len(n), each = 4L)
    x1    <- rnorm(length(t_vec))
    y     <- round(0.5 * x1 + rnorm(length(t_vec), sd = 0.4), 1L)
    df    <- data.frame(y = y, x1 = x1, subj = subj, t = t_vec)

    # initial = c(intercept_coef, x1_coef) as a unit direction
    fit_SensIAT_single_index_fixed_coef_model(
        y ~ x1,
        data    = df,
        id      = subj,
        kernel  = "K2_Biweight",
        initial = c(0, 1)          # unit vector along x1 direction
    )
}

test_that("model object caches Xb_train and Yi at fit time", {
    m <- make_test_model()
    expect_false(is.null(m$Xb_train), label = "Xb_train is cached")
    expect_false(is.null(m$Yi),       label = "Yi is cached")

    # Manually recompute and compare
    Xi_check    <- model.matrix(terms(m), m$data)
    Xb_expected <- as.vector(Xi_check %*% m$coef)
    Yi_expected <- model.response(model.frame(m))

    expect_equal(m$Xb_train, Xb_expected)
    expect_equal(m$Yi,       Yi_expected)
})

test_that("compute_SensIAT_expected_values matches reference for single query, single alpha", {
    m     <- make_test_model()
    nd    <- m$data[1L, , drop = FALSE]
    alpha <- 0.1

    result <- compute_SensIAT_expected_values(m, alpha = alpha, new.data = nd)

    # Reference via pcoriaccel_NW CDF path
    Xi_new <- model.matrix(terms(m), data = nd)
    xb_q   <- as.vector(Xi_new %*% m$coef)
    ref    <- ref_expectations(m$Xb_train, m$Yi, xb_q, alpha, m$bandwidth,
                                kernel = attr(m, "kernel"))

    expect_equal(result$E_exp_alphaY,  ref$E_exp,  tolerance = 1e-10)
    expect_equal(result$E_Yexp_alphaY, ref$E_Yexp, tolerance = 1e-10)
})

test_that("compute_SensIAT_expected_values matches reference for multiple query rows", {
    m     <- make_test_model()
    nd    <- m$data[1:5, , drop = FALSE]
    alpha <- -0.05

    result <- compute_SensIAT_expected_values(m, alpha = alpha, new.data = nd)

    Xi_new <- model.matrix(terms(m), data = nd)
    xb_q   <- as.vector(Xi_new %*% m$coef)
    ref    <- ref_expectations(m$Xb_train, m$Yi, xb_q, alpha, m$bandwidth,
                                kernel = attr(m, "kernel"))

    expect_equal(result$E_exp_alphaY,  ref$E_exp,  tolerance = 1e-10)
    expect_equal(result$E_Yexp_alphaY, ref$E_Yexp, tolerance = 1e-10)
})

test_that("compute_SensIAT_expected_values handles alpha = 0 (moments reduce to 1 and E[Y])", {
    m     <- make_test_model()
    nd    <- m$data[1:3, , drop = FALSE]

    result <- compute_SensIAT_expected_values(m, alpha = 0, new.data = nd)

    # When alpha=0: exp(0*Y)=1 for all Y, so E[exp(0*Y)]=1, E[Y*exp(0*Y)]=E[Y]
    expect_equal(result$E_exp_alphaY, rep(1, 3), tolerance = 1e-10)
    # E_Y_past = E_Yexp / E_exp = E[Y]
    expect_equal(result$E_Y_past, result$E_Yexp_alphaY, tolerance = 1e-12)
})

test_that("compute_SensIAT_expected_values handles vector alpha (dispatch path)", {
    m      <- make_test_model()
    nd     <- m$data[1:2, , drop = FALSE]
    alphas <- c(-0.1, 0, 0.1)

    result <- compute_SensIAT_expected_values(m, alpha = alphas, new.data = nd)

    expect_equal(nrow(result), 2L * length(alphas))
    for (a in alphas) {
        ref_a <- compute_SensIAT_expected_values(m, alpha = a, new.data = nd)
        rows  <- result[result$alpha == a, , drop = FALSE]
        expect_equal(rows$E_exp_alphaY,  ref_a$E_exp_alphaY,  tolerance = 1e-12)
        expect_equal(rows$E_Yexp_alphaY, ref_a$E_Yexp_alphaY, tolerance = 1e-12)
    }
})

test_that("pcoriaccel_NW_expectations matches ref_expectations directly", {
    set.seed(42L)
    n     <- 200L
    Xb    <- rnorm(n)
    Yi    <- round(rnorm(n, sd = 0.5), 1L)
    xb_q  <- seq(-1.5, 1.5, length.out = 20L)
    alpha <- 0.2
    h     <- 0.5

    mat <- pcoriaccel_NW_expectations(Xb, Yi, xb_q, alpha, h, "K2_Biweight")
    ref <- ref_expectations(Xb, Yi, xb_q, alpha, h, "K2_Biweight")

    expect_equal(mat[, 1L], ref$E_exp,  tolerance = 1e-10)
    expect_equal(mat[, 2L], ref$E_Yexp, tolerance = 1e-10)
})

test_that("pcoriaccel_NW_expectations matches for dnorm kernel", {
    set.seed(7L)
    n     <- 150L
    Xb    <- rnorm(n)
    Yi    <- round(rnorm(n, sd = 0.4), 1L)
    xb_q  <- seq(-1, 1, length.out = 10L)
    alpha <- -0.1
    h     <- 0.8

    mat <- pcoriaccel_NW_expectations(Xb, Yi, xb_q, alpha, h, "dnorm")
    ref <- ref_expectations(Xb, Yi, xb_q, alpha, h, "dnorm")

    expect_equal(mat[, 1L], ref$E_exp,  tolerance = 1e-10)
    expect_equal(mat[, 2L], ref$E_Yexp, tolerance = 1e-10)
})
