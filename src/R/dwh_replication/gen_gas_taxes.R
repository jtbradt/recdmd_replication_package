############################## FRONT MATTER ##############################
#   SCRIPT:     gen_gas_taxes.R
#   AUTHOR:     Jacob Bradt (jbradt@g.harvard.edu)
#   NOTES:      Script to construct destination state gas tax data

# Load dependencies:
pacman::p_load(here, data.table, tidyverse, fixest, sf, units, zoo, fredr, tigris, haven)

# Set project root:
i_am("R/gen_gas_taxes.R")

# STEP 1: Import and format sites data ----

# Import trips data to extract destination lat-long:
trips_local <- fread(here("data/trips_local.csv"))
sites <- distinct(trips_local[, c("site_id_agg", "lat_lon_des_agg")])
sites <- sites[site_id_agg != -1]
rm(trips_local)

# Format site lat-long:
sites[, (c("lat", "long")) := list(as.numeric(substr(lat_lon_des_agg, 1, 7)), -as.numeric(substr(lat_lon_des_agg, 10, 16)))]
sites[, lat_lon_des_agg := NULL]

# Create site spatial data:
sites_sf <- st_as_sf(sites, coords = c("long", "lat"),crs="EPSG:4326") %>%
  st_transform(., st_crs("ESRI:102005"))

# STEP 2: Import and format respondent data ----

# Import and format local respondents data:
respondents_local <- fread(here("data/respondents_local.csv"))
respondents_local <- distinct(respondents_local[, c("caseid", "state")])

# Import and format national respondents data:
respondents_national <- fread(here("data/respondents_national.csv"))
respondents_national <- distinct(respondents_national[, c("caseid", "state")])

# Merge national and local origins data:
respondents <- rbindlist(list(respondents_local, respondents_national), fill = T)
rm(respondents_local, respondents_national)

# Convert state abbr to fips code:
respondents <- merge(respondents, distinct(fips_codes[, 1:2]), by = "state")
respondents[, state_code := as.numeric(state_code)]
respondents[, state := NULL]

# STEP 3: Assign sites to states ----

# Import states shapefile:
states_sf <- tigris::states() %>%
  st_as_sf(.) %>%
  st_transform(., st_crs("ESRI:102005"))

# Construct state distance matrix:
states_distmat <- st_distance(st_centroid(states_sf), st_centroid(states_sf))
rownames(states_distmat) <- states_sf$STATEFP
colnames(states_distmat) <- states_sf$STATEFP

# Assign sites to states
sites_sf <- st_join(sites_sf, states_sf[, c("STATEFP")])
sites <- st_drop_geometry(sites_sf)
rm(states_sf, sites_sf)

# Format sites data:
sites <- setDT(sites)[order(site_id_agg)]
sites[, state := as.numeric(STATEFP)]
sites[, STATEFP := NULL]

# STEP 4: Import state gas tax data ----

# Import state gas tax data:
state_taxes <- read_dta(here("data/annual_taxes_cleaned.dta")) %>%
  data.table(.)

# Subset to 2012 and merge to sites:
state_taxes <- state_taxes[ year == 2012]
sites <- merge(sites, state_taxes[, c("fips", "gas_tax")], by.x = "state", by.y = "fips")
setnames(sites, "gas_tax", "dest_gas_tax")

# Subset to 2012 and merge to origins:
respondents <- merge(respondents, state_taxes[, c("fips", "gas_tax")], by.x = "state_code", by.y = "fips")
setnames(respondents, "gas_tax", "origin_gas_tax")

# Expand origins data to long:
respondents <- merge(respondents, CJ(caseid = respondents$caseid, site_id_agg = 1:83), by = "caseid")

# Merge destination gas tax data to origins long:
respondents <- merge(respondents, sites, by = "site_id_agg")
rm(sites)

# STEP 5: Construct weighted measure of intermediate states gas tax ----

# Ordered state gas taxes:
state_taxes_order <- lapply(colnames(states_distmat), function(x){
  if (as.numeric(x) %in% state_taxes$fips) {
    state_taxes[fips == as.numeric(x)]$gas_tax
  } else{
    0.0
  }
}) %>% do.call(rbind, .)

# Calculated average other state gas taxes weighted by inverse distance:
states_distmat = 1 / (states_distmat/1000)
diag(states_distmat) <- 0.0
state_taxes_int <- data.table(state_code = as.numeric(rownames(states_distmat)), int_gas_tax = (states_distmat %*% (state_taxes_order))[,1])
rm(state_taxes_order, states_distmat)

# Merge intermediate gas tax data to full data:
respondents <- merge(respondents, state_taxes_int, by = "state_code")
rm(state_taxes_int, state_taxes)

# Save destination state gas tax data:
respondents[, (c("state_code", "state")) := NULL]
fwrite(respondents, here("data/gas_tax_delta.csv"))
rm(respondents)
