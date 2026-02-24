# Debug time_var type
library(dplyr)
devtools::load_all("/workspaces/pcoriRPackage")

# Simulate what the main function does
data <- SensIAT_example_data %>%
  group_by(Subject_ID) %>%
  mutate(
    ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - lag(Time, default = NA_real_, order_by = Time)
  ) %>%
  ungroup() %>%
  filter(Subject_ID <= 5)

# Simulate enquo(time) as it happens in main function
time <- data$Time
time.var <- rlang::enquo(time)

cat("time.var class:", class(time.var), "\n")
cat("time.var is_quosure:", rlang::is_quosure(time.var), "\n")
cat("quo_name:", rlang::quo_name(time.var), "\n")

# Try with the actual call
test_fn <- function(time) {
  time.var <- rlang::enquo(time)
  cat("Inside function - time.var class:", class(time.var), "\n")
  cat("Inside function - is_quosure:", rlang::is_quosure(time.var), "\n")
  cat("Inside function - quo_name:", rlang::quo_name(time.var), "\n")
  time.var
}

tv <- test_fn(data$Time)
cat("\nReturned time.var class:", class(tv), "\n")
cat("Returned is_quosure:", rlang::is_quosure(tv), "\n")
cat("quo_name on returned:", rlang::quo_name(tv), "\n")

# What does eval_tidy do with this?
cat("\neval_tidy(tv, data):", head(rlang::eval_tidy(tv, data)), "\n")
