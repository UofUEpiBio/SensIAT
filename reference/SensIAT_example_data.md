# SensIAT Example Data

A simulated dataset for use in the SensIAT tutorial, testing and
documentation.

## Usage

``` r
SensIAT_example_data

SensIAT_example_fulldata
```

## Format

`SensIAT_example_data` is a data frame with 779 rows and 4 variables
consisting of 200 simulated patients. Each row in the data represents a
visit for the patient. The columns are:

- Subject_ID:

  A unique identifier for each patient.

- Visit:

  The ordinal number of the visit for the patient. Baseline observation
  is 0.

- Time:

  The time of the visit in days, since baseline.

- Outcome:

  The outcome of interest.

`SensIAT_example_fulldata` is a data frame with 1614 rows and 5
variables consisting of 400 simulated patients, 200 for each treatment
arm. Each row in the data represents a visit for the patient. The
columns are:

- Subject_ID:

  A unique identifier for each patient.

- Visit:

  The ordinal number of the visit for the patient. Baseline observation
  is 0.

- Time:

  The time of the visit in days, since baseline.

- Outcome:

  The outcome of interest.

- Treatment_group:

  Treatment or control group.

## Functions

- `SensIAT_example_fulldata`: A simulated dataset with both treatment
  and control groups.
