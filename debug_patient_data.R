# Debug patient_data issue
library(dplyr)
devtools::load_all("/workspaces/pcoriRPackage")

set.seed(456)
data_with_lags <- SensIAT_example_data %>%
  group_by(Subject_ID) %>%
  mutate(
    ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - lag(Time, default = NA_real_, order_by = Time)
  ) %>%
  ungroup() %>%
  filter(Subject_ID <= 5)

# Check what data looks like
cat("=== Original data structure ===\n")
print(head(data_with_lags))

# What main function would see
data <- data_with_lags
id <- data$Subject_ID
unique_ids <- unique(id)

# How main function builds patient_data_list
patient_data_list <- lapply(unique_ids, function(pid) {
    data[id == pid, ]
})
names(patient_data_list) <- as.character(unique_ids)

cat("\n=== Patient 1 from patient_data_list ===\n")
print(patient_data_list[["1"]])

# How my debug script builds patient data
cat("\n=== Patient 1 from my debug script ===\n")
patient_data_manual <- data_with_lags %>% filter(Subject_ID == 1)
print(patient_data_manual)

# Compare time ranges
cat("\n=== Time ranges ===\n")
cat("From patient_data_list:", range(patient_data_list[["1"]]$Time), "\n")
cat("From manual filter:", range(patient_data_manual$Time), "\n")

# Check column names
cat("\n=== Column names ===\n")
cat("From patient_data_list:", names(patient_data_list[["1"]]), "\n")
cat("From manual filter:", names(patient_data_manual), "\n")
