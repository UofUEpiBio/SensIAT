test_that("add_class adds S3 class to object", {
    # Test adding class to a list
    obj <- list(a = 1, b = 2)
    result <- add_class(obj, "my_class")
    
    expect_s3_class(result, "my_class")
    expect_s3_class(result, "list")
    expect_equal(class(result), c("my_class", "list"))
})

test_that("add_class preserves existing classes", {
    # Create object with existing class
    obj <- structure(list(x = 1), class = c("existing", "list"))
    result <- add_class(obj, "new_class")
    
    expect_equal(class(result), c("new_class", "existing", "list"))
})

test_that("add_class doesn't duplicate classes", {
    # Add class that already exists
    obj <- structure(list(x = 1), class = c("my_class", "list"))
    result <- add_class(obj, "my_class")
    
    expect_equal(class(result), c("my_class", "list"))
    # Should not have duplicated "my_class"
    expect_equal(sum(class(result) == "my_class"), 1)
})

test_that("add_class works with different object types", {
    # Test with data.frame
    df <- data.frame(x = 1:3)
    result_df <- add_class(df, "custom_df")
    expect_s3_class(result_df, "custom_df")
    expect_s3_class(result_df, "data.frame")
    
    # Test with vector
    vec <- c(1, 2, 3)
    result_vec <- add_class(vec, "custom_vector")
    expect_s3_class(result_vec, "custom_vector")
    
    # Test with logical
    lgl <- TRUE
    result_lgl <- add_class(lgl, "flag")
    expect_s3_class(result_lgl, "flag")
})

test_that("add_class validates input", {
    obj <- list(a = 1)
    
    # Should error if class is not character
    expect_error(
        add_class(obj, 123),
        "class must be a character vector"
    )
    
    expect_error(
        add_class(obj, NULL),
        "class must be a character vector"
    )
})

test_that("add_class works with S3 objects only", {
    # This test documents that add_class is for S3 only
    # S3 classes work with simple class() assignment
    obj <- list(value = 42)
    result <- add_class(obj, "S3_class")
    
    expect_true(inherits(result, "S3_class"))
    expect_equal(class(result)[1], "S3_class")
    
    # Note: This function does NOT work with S4 or R6 classes
    # S4 would require setClass/new, R6 would require R6Class
})

test_that("add_class maintains object attributes", {
    # Create object with attributes
    obj <- list(x = 1)
    attr(obj, "custom_attr") <- "test_value"
    
    result <- add_class(obj, "new_class")
    
    expect_equal(attr(result, "custom_attr"), "test_value")
    expect_equal(result$x, 1)
})
