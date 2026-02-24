# Test different ways to extract column name from quosure
library(rlang)

# Case 1: User passes unquoted column name
test1 <- function(time) {
  time.var <- enquo(time)
  cat("Case 1 - quo_get_expr:", deparse(quo_get_expr(time.var)), "\n")
}

# Case 2: User passes df$Time  
test2 <- function(time) {
  time.var <- enquo(time)
  cat("Case 2 - quo_get_expr:", deparse(quo_get_expr(time.var)), "\n")
}

# Simulate data
df <- data.frame(Time = 1:10, X = 11:20)

cat("=== Testing different call styles ===\n\n")

cat("Calling with Time (symbol):\n")
test1(Time)

cat("\nCalling with df$Time (expression):\n")
test2(df$Time)

# How to extract column name from either case
extract_colname <- function(time.var) {
  expr <- rlang::quo_get_expr(time.var)
  
  # Case 1: It's a symbol (unquoted column name)
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  
  # Case 2: It's a call like df$Time
  if (is.call(expr) && identical(expr[[1]], as.symbol("$"))) {
    return(as.character(expr[[3]]))
  }
  
  # Case 3: It's a literal vector - try to find matching column
  return(NULL)
}

cat("\n=== Testing extract_colname ===\n")

test_extract <- function(time) {
  time.var <- enquo(time)
  result <- extract_colname(time.var)
  cat("Extracted column name:", result, "\n")
  result
}

cat("With Time:\n")
test_extract(Time)

cat("\nWith df$Time:\n")
test_extract(df$Time)
