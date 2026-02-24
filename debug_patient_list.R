# Debug: check patient_data_list contents
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

# Simulate main function's patient_data_list creation
data <- data_with_lags
id <- data$Subject_ID
unique_ids <- unique(id)

patient_data_list <- lapply(unique_ids, function(pid) {
    data[id == pid, ]
})
names(patient_data_list) <- as.character(unique_ids)

cat("Patient data list contents:\n")
for (key in names(patient_data_list)) {
  pd <- patient_data_list[[key]]
  cat("Key:", key, "nrow:", nrow(pd), "Times:", paste(pd$Time, collapse=", "), "\n")
}
