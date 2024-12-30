############################## FRONT MATTER ##############################
#   SCRIPT:     gen_census_vars.R
#   AUTHOR:     Jacob Bradt (jbradt@g.harvard.edu)
#   NOTES:      

# Load dependencies:
pacman::p_load(here, data.table, tidyverse, tidycensus)

# Set project root:
i_am("R/gen_census_vars.R")

# STEP 1: Import ACS 5-year data ----

# Get ACS data from Census API:
acs_2012 <- get_acs("cbsa", variables = c("B01003_001", "B06011_001", "B23025_002", "B25103_001", "B23025_004"), year = 2012) %>%
  data.table(.)

# Convert from long to wide:
acs_2012 <- dcast(acs_2012, GEOID  ~ variable, value.var = "estimate")

# Split ZCTA and state fields:
#acs_2012[, (c("GEOID", "state")) := list(substr(GEOID, 3, 7), substr(GEOID, 1, 2))]

# Set field names:
setnames(
  acs_2012,
  c(
    "GEOID",
    "B01003_001",
    "B06011_001",
    "B23025_002",
    "B23025_004",
    "B25103_001"
  ),
  c("cbsa", "pop", "med_income", "laborforce", "employed", "proptax")
)

# Construct employment rate:
acs_2012[, employment_rate := employed / laborforce]


# STEP 2: Save ACS data ----

fwrite(acs_2012, here("data/acs5yr_2012.csv"))
rm(acs_2012)

