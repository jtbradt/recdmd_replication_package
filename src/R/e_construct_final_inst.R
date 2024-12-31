############################## FRONT MATTER ##############################
#   SCRIPT:     construct_final_inst.R
#   AUTHOR:     Jacob Bradt (jacob.bradt@mccombs.utexas.edu)
#   NOTES:      Script to construct first stage cost instruments for
#               replicating English et al (2018) Deepwater Horizon study.
#               Also estimates versions of first stage travel cost regressions
#               reported in paper.

# Load dependencies:
pacman::p_load(here, data.table, tidyverse, fixest, sf, units, zoo, fredr, tidycensus, haven)

# Set project root:
i_am("src/R/e_construct_final_inst.R")

# STEP 1: Import distance matrix ----

# Import distance matrix:
origins <- fread(here("data/dwh_replication/distance_mat.csv"))

# STEP 2: Import and merge travel cost data ----

# Travel cost data:
tc <- fread(here("data/dwh_replication/tc1.txt"))

# Setnames:
names(tc) <- c("caseid", "es_origin", "income", paste0("tc_", 1:83))

# Establish caseid order:
caseids <- tc[, 1]

# Merge travel costs to origins data:
origins <- merge(origins, tc, by = c("caseid", "es_origin"))
rm(tc)

# Import choices:
choices <- fread(here("data/dwh_replication/trips_demos.txt"))
choices <- choices[, 1:3]
names(choices) <- c("caseid", "es_origin", "co")
summary(choices$caseid == caseids$caseid)

# Identify local observations:
local_respondents <- fread(here("data/dwh_replication/respondents_local.csv"))
choices[, local := ifelse(caseid %in% local_respondents$caseid, 1, 0)]
rm(local_respondents)

# Construct weight and merge choice occasion data to full origins data:
choices[, wgt := co / sum(co)]
choices[, (c("wgt_local", "wgt_ntl")) := list(ifelse(local == 1, co / sum(choices[local == 1]$co), 0),
                                              ifelse(local == 0, co / sum(choices[local == 0]$co), 0))]
origins <- merge(origins, choices)
rm(choices)
summary(origins$caseid == caseids$caseid)

# STEP 3: Import and merge WTI data ----

# Import WTI data:
wti <- fread(here("data/dwh_replication/eia_wti.csv"))
wti[,(c("year", "mon")) := list(year(as.yearmon(date, "%b-%Y")), month(as.yearmon(date, "%b-%Y")))]
wti[, date := NULL]

# Convert WTI data to 2011 dollars:
fredr::fredr_set_key("588eb023c04a4e8fd31b6b828b6ddbd6")
cpi <- fredr::fredr("CPIAUCSL") %>%
  data.table(.)
cpi[, (c("year", "mon")) := list(year(date), month(date))]
cpi <- cpi[, .(cpi = mean(value)), by = c("year", "mon")]
wti <- merge(wti, cpi)
cpi_2011 <- mean(wti[year == 2011]$cpi)
wti[, wti := wti * cpi_2011 / cpi]
wti[, cpi := NULL]
rm(cpi, cpi_2011)

# Merge WTI to origins data:
origins <- merge(origins, wti, by = c("year", "mon"), sort=FALSE)
summary(origins$caseid == caseids$caseid)
rm(wti)

# Construct wide instrument matrix:
Zmat1 <- origins[, grepl("dist",names(origins)), with = F] * origins$wti
write.table(
  Zmat1,
  here("data/dwh_replication/zmat/zmat01.txt"),
  row.names = F,
  col.names = F
)
rm(Zmat1)

# STEP 4: Import and merge oil supply shock data ----

# Import estimated structural shocks to oil supply:
supply_shocks <- fread(here("data/dwh_replication/oil_struc_shocks.csv"))
supply_shocks[, year_mon := as.yearmon(year_mon)]
supply_shocks[, (c("year", "mon")) := list(year(year_mon), month(year_mon))]

# Merge structural oil supply shocks to origins data:
origins <- merge(origins, supply_shocks, by = c("year", "mon"), sort=FALSE)
summary(origins$caseid == caseids$caseid)
rm(supply_shocks)

# Construct wide instrument matrix:
Zmat2 <- origins[, grepl("dist",names(origins)), with = F] * origins$d_prod
write.table(
  Zmat2,
  here("data/dwh_replication/zmat/zmat02.txt"),
  row.names = F,
  col.names = F
)
Zmat3 <- origins[, grepl("dist",names(origins)), with = F] * origins$refiner_price
write.table(
  Zmat3,
  here("data/dwh_replication/zmat/zmat03.txt"),
  row.names = F,
  col.names = F
)
rm(Zmat2, Zmat3)

# STEP 5: Import and merge commuting time data ----

# Import respondents avg. commute time:
respondents_commute <- fread(here("data/dwh_replication/respondents_commutes.csv"))

# Add CBSA:
cbsa <- fread(here("data/dwh_replication/ZIP_CBSA_092012.csv"))
respondents_commute <- merge(respondents_commute, cbsa[, 1:3], by.x = "zipcode", by.y = "ZIP")
respondents_commute <- respondents_commute[respondents_commute[, .I[which.max(RES_RATIO)], by=c("zipcode", "caseid", "cutoffym")]$V1]
rm(cbsa)

# Merge commute times to origins data:
origins <- merge(origins, respondents_commute[, c(1, 2, 4, 5)], by = "caseid", all.x = T, sort = F)
origins[, mean_commute := ifelse(is.na(mean_commute), mean(origins$mean_commute, na.rm = T), mean_commute)]
origins[, CBSA := ifelse(is.na(CBSA), 99999, CBSA)]
rm(respondents_commute)
summary(origins$caseid == caseids$caseid)

# Construct wide instrument matrix:
Zmat4 <- matrix(rep(origins$mean_commute, sum(grepl("dist",names(origins)))), ncol = sum(grepl("dist",names(origins))))
write.table(
  Zmat4,
  here("data/dwh_replication/zmat/zmat04.txt"),
  row.names = F,
  col.names = F
)
rm(Zmat4)

# STEP 6: Import and merge gas tax data ----

# Reshape to long format:
long_vars <- list(names(origins)[grepl("dist_", names(origins))], names(origins)[grepl("tc_", names(origins))])
origins_long <- melt(origins, measure.vars = long_vars, value.name = c("dist", "tc"))
rm(long_vars)

# Import destination tax data:
gas_tax <- fread(here("data/dwh_replication/gas_tax_delta.csv"))

# Reshape destination gas tax to wide:
gas_tax_wide <- dcast(gas_tax, caseid ~ site_id_agg, value.var = "dest_gas_tax")
gas_tax_wide <- merge(caseids, gas_tax_wide, by = "caseid", all.x = T)
summary(gas_tax_wide$caseid == caseids$caseid)

# Construct wide instrument matrix for destination gas tax:
Zmat5 <- origins[, grepl("dist",names(origins)), with = F] * gas_tax_wide[, 2:ncol(gas_tax_wide)]
write.table(
  Zmat5,
  here("data/dwh_replication/zmat/zmat05.txt"),
  row.names = F,
  col.names = F
)
rm(Zmat5, gas_tax_wide)

# Reshape intermediate gas tax to wide:
gas_tax_wide <- dcast(gas_tax, caseid ~ site_id_agg, value.var = "int_gas_tax")
gas_tax_wide <- merge(caseids, gas_tax_wide, by = "caseid", all.x = T)
summary(gas_tax_wide$caseid == caseids$caseid)

# Construct wide instrument matrix for intermediate gas tax:
Zmat6 <- origins[, grepl("dist",names(origins)), with = F] * gas_tax_wide[, 2:ncol(gas_tax_wide)]
write.table(
  Zmat6,
  here("data/dwh_replication/zmat/zmat06.txt"),
  row.names = F,
  col.names = F
)
rm(Zmat6, gas_tax_wide)

# Merge to long data:
origins_long[, variable := as.numeric(variable)]
origins_long <- merge(origins_long, gas_tax, by.x = c("caseid", "variable"), by.y = c("caseid", "site_id_agg"), all.x = T)
rm(gas_tax)

# STEP 7: CBSA field for use in FS estimation ----

# Add state field:
respondents_local <- fread(here("data/dwh_replication/respondents_local.csv"))
respondents_national <- fread(here("data/dwh_replication/respondents_national.csv"))
respondents <- rbind(distinct(respondents_local[, c("caseid", "state", "edu_cat")]), distinct(respondents_national[, c("caseid", "state", "edu_cat")]))
origins <- merge(origins, respondents)
origins_long <- merge(origins_long, respondents, by = "caseid")
rm(respondents_local, respondents_national, respondents)

# Extract CBSA field:
summary(origins$caseid == caseids$caseid)
write.table(
  origins$state,
  here("data/dwh_replication/zmat/zmat07.txt"),
  row.names = F,
  col.names = F
)

# STEP 8: Merge ACS 5-year data ----

# Import ACS 5-yr estimates at ZCTA level for 2012:
acs_2012 <- fread(here("data/dwh_replication/acs5yr_2012.csv"))

# Merge ACS data to wide origins data:
origins <- merge(origins, acs_2012, by.x = "CBSA", by.y = "cbsa", all.x = T)
origins <- origins[order(caseid)]
summary(origins$caseid == caseids$caseid)

# Remove NAs from ACS data:
acs_vars <- c("pop", "med_income", "laborforce", "proptax", "employment_rate")
origins[, (acs_vars) := lapply(.SD, function(x){ifelse(is.na(x), mean(x, na.rm = T), x)}), .SDcols = acs_vars]

# Save ACS data as instrument matrix:
Zmat8 <- matrix(rep(log(origins$med_income), sum(grepl("dist",names(origins)))), ncol = sum(grepl("dist",names(origins))))
write.table(
  Zmat8,
  here("data/dwh_replication/zmat/zmat08.txt"),
  row.names = F,
  col.names = F
)
Zmat9 <-matrix(rep(log(origins$employment_rate), sum(grepl("dist",names(origins)))), ncol = sum(grepl("dist",names(origins))))
write.table(
  Zmat9,
  here("data/dwh_replication/zmat/zmat09.txt"),
  row.names = F,
  col.names = F
)
Zmat10 <- matrix(rep(log(origins$pop), sum(grepl("dist",names(origins)))), ncol = sum(grepl("dist",names(origins))))
write.table(
  Zmat10,
  here("data/dwh_replication/zmat/zmat10.txt"),
  row.names = F,
  col.names = F
)
Zmat11 <- as.matrix(origins[, grepl("dist",names(origins)), with = F])
write.table(
  Zmat11,
  here("data/dwh_replication/zmat/zmat11.txt"),
  row.names = F,
  col.names = F
)
rm(Zmat8, Zmat9, Zmat10, Zmat11, origins, caseids)

# Merge to long origins data:
origins_long <- merge(origins_long, acs_2012, by.x = "CBSA", by.y = "cbsa", all.x = T)

# Assign missing values:
origins_long[, (acs_vars) := lapply(.SD, function(x){ifelse(is.na(x), mean(x, na.rm = T), x)}), .SDcols = acs_vars]
rm(acs_2012, acs_vars)

# STEP 8: Estimate first stage for reporting in text ----

# First stage estimation:
origins_long[, employment_rate := employment_rate * 100]
setFixest_fml(..fs_fml = ~ (dest_gas_tax + int_gas_tax) * dist - dist - dest_gas_tax - int_gas_tax)
fs <- list(
  #feols(tc ~ wti * dist - wti + ..fs_fml | variable, data = origins_long),
  #feols(tc ~ (d_prod + refiner_price) * dist - refiner_price - d_prod  + ..fs_fml | variable, data = origins_long),
  #feols(tc ~ (d_prod + refiner_price + log(employment_rate) + log(med_income) + log(pop)) * dist - d_prod - refiner_price + ..fs_fml | variable, data = origins_long),
  #feols(tc ~ wti * dist - wti + ..fs_fml | variable, weights = ~wgt, data = origins_long),
  feols(tc ~ (d_prod + refiner_price) * dist - refiner_price - d_prod  + ..fs_fml | variable, weights = ~wgt, data = origins_long),
  feols(tc ~ (d_prod + refiner_price) * dist + log(employment_rate) + log(med_income) + log(pop) - d_prod - refiner_price + ..fs_fml | variable, weights = ~wgt, data = origins_long),
  feols(tc ~ (d_prod + refiner_price + log(employment_rate) + log(med_income) + log(pop)) * dist - d_prod - refiner_price + ..fs_fml | variable, weights = ~wgt, data = origins_long)
)

# Local vs. national first stage estimation:
fs_local <- list(
  feols(tc ~ (d_prod + refiner_price + log(employment_rate) + log(med_income) + log(pop)) * dist - d_prod - refiner_price + ..fs_fml | variable, weights = ~wgt, data = origins_long),
  feols(tc ~ (d_prod + refiner_price + log(employment_rate) + log(med_income) + log(pop)) * dist - d_prod - refiner_price + ..fs_fml | variable, weights = ~wgt_ntl, data = origins_long[local == 0]),
  feols(tc ~ (d_prod + refiner_price + log(employment_rate) + log(med_income) + log(pop)) * dist - d_prod - refiner_price + ..fs_fml | variable, weights = ~wgt_local, data = origins_long[local == 1])
)

# Save first stage estimates:
setFixest_dict(
  c(
    "tc" = "Travel Cost",
    "wti" = "WTI",
    "dist" = "Distance",
    "gas_tax_delta" = "Gas Tax Diff.",
    "dest_gas_tax" = "Destination Gas Tax",
    "int_gas_tax" = "Mean Route Gas Tax",
    "mean_commute" = "Origin Commute Time",
    "d_prod" = "Oil Supply Shock",
    "refiner_price" = "Oil Demand Shock",
    "variable" = "Alternative Specific Constants",
    "log(emloyment_rate)" = "log(Employment Rate)",
    "log(med_income)" = "log(Median Income)",
    "log(pop)" = "log(Population)"
  )
)
etable(
  fs,
  digits = 3,
  tex = T,
  powerBelow = -6,
  se.below = F,
  digits.stats = 3,
  signif.code = NA,
  fitstat = ~ n + r2 + wr2 + wald,
  extralines = list("_^Sample Weights" = c("No", "No", "Yes", "Yes")),
  file = here("output/dwh_replication/dwh_first_stage_main.tex"),
  replace = T
)
etable(
  fs_local,
  digits = 3,
  tex = T,
  powerBelow = -6,
  se.below = F,
  digits.stats = 3,
  signif.code = NA,
  fitstat = ~ n + r2 + wr2 + wald,
  extralines = list("_^Sample Weights" = c("Yes", "Yes", "Yes"),
                    "_^Sample" = c("National", "Non-local Respondents", "Local Respondents")),
  file = here("output/dwh_replication/dwh_first_stage_local.tex"),
  replace = T
)
rm(fs, fs_local, origins_long)
