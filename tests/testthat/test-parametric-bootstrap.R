library(testthat)

test_that("sample_parametric_coeffs works on lm objects", {
  set.seed(1)
  df <- data.frame(x = rnorm(20))
  df$y <- 1 + 2 * df$x + rnorm(20, 0, 0.1)
  m <- lm(y ~ x, data = df)
  samp <- SensIAT:::sample_parametric_coeffs(m)
  expect_true(is.numeric(samp))
  expect_equal(length(samp), length(coef(m)))
  expect_true(all(names(samp) %in% names(coef(m))))
})

test_that("parametric_bootstrap runs with a simple coxph intensity model", {
  skip_if_not_installed("survival")
  # create simple survival data to fit a coxph
  set.seed(2)
  n <- 50
  df <- data.frame(x = rnorm(n))
  df$time <- rexp(n, rate = exp(-0.5 * df$x))
  df$status <- rbinom(n, 1, 0.9)
  fit <- survival::coxph(survival::Surv(time, status) ~ x, data = df)

  res <- SensIAT:::parametric_bootstrap(
    nboot = 2,
    intensity_model = fit,
    outcome_model = NULL,
    simulate_args = list(
      n_subjects = 3,
      End = 5,
      max_visits = 5
    ),
    seed = 42
  )

  expect_length(res, 2)
  expect_true(all(vapply(res, function(x) inherits(x, "data.frame"), logical(1))))
})

test_that("parametric_bootstrap_within_group runs on a small within-group model", {
  skip_if_not_installed("survival")
  data("SensIAT_example_data", package = "SensIAT", envir = environment())
  fixture_data <- dplyr::filter(SensIAT_example_data, Subject_ID %in% head(unique(SensIAT_example_data$Subject_ID), 8))

  model <- fit_SensIAT_within_group_model(
    group.data = fixture_data,
    outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
    alpha = 0,
    id = Subject_ID,
    outcome = Outcome,
    time = Time,
    End = 830,
    knots = c(60, 260, 460)
  )

  res <- SensIAT:::parametric_bootstrap_within_group(
    nboot = 1,
    within_group_model = model,
    simulate_args = list(
      n_subjects = 3,
      End = 5,
      max_visits = 5
    ),
    seed = 123
  )

  expect_length(res, 1)
  expect_true(inherits(res[[1]], "data.frame"))
})

test_that("single-index outcome models support vcov() and coefficient sampling", {
  data("SensIAT_example_data", package = "SensIAT", envir = environment())
  fixture_data <- dplyr::filter(SensIAT_example_data, Subject_ID %in% head(unique(SensIAT_example_data$Subject_ID), 8))

  model <- fit_SensIAT_within_group_model(
    group.data = fixture_data,
    outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
    alpha = 0,
    id = Subject_ID,
    outcome = Outcome,
    time = Time,
    End = 830,
    knots = c(60, 260, 460)
  )

  outcome_model <- model$models$outcome
  vc <- vcov(outcome_model)

  expect_true(is.matrix(vc))
  expect_equal(nrow(vc), length(coef(outcome_model)))
  expect_equal(ncol(vc), length(coef(outcome_model)))

  samp <- SensIAT:::sample_parametric_coeffs(outcome_model)
  expect_length(samp, length(coef(outcome_model)))
  expect_equal(names(samp), names(coef(outcome_model)))
})
