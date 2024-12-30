############################## FRONT MATTER ##############################
#   SCRIPT:     gen_dist_mat.R
#   AUTHOR:     Jacob Bradt (jbradt@g.harvard.edu)
#   NOTES:      Script to construct origin-desination distance matrix for
#   replicating English et al (2018) Deepwater Horizon study

# Load dependencies:
pacman::p_load(here, data.table, tidyverse, fixest, sf, units, zoo, fredr, tidycensus, tigris, ggspatial, cowplot)

# Set project root:
i_am("R/gen_instrument.R")

# STEP 1: Import and format sites data ----

# Import trips data to extract destination lat-long:
trips_local <- fread(here("data/trips_local.csv"))
sites <- distinct(trips_local[, c("site_id_agg", "lat_lon_des_agg")])
sites <- sites[site_id_agg != -1]
rm(trips_local)

# Format site lat-long:
sites[, (c("lat", "long")) := list(as.numeric(substr(lat_lon_des_agg, 1, 7)), -as.numeric(substr(lat_lon_des_agg, 10, 16)))]
sites[, lat_lon_des_agg := NULL]

# STEP 2: Import and format respondent data ----

# Import and format local respondents data:
respondents_local <- fread(here("data/respondents_local.csv"))
respondents_local <- distinct(respondents_local[, c("caseid", "cutoffym")])
origins_local <- fread(here("data/origins_local.csv"))
origins_local <- distinct(origins_local[, c("caseid", "es_origin", "lat_lon_origin_final")])
origins_local <- merge(origins_local, respondents_local)
rm(respondents_local)

# Import and format national respondents data:
respondents_national <- fread(here("data/respondents_national.csv"))
respondents_national <- distinct(respondents_national[, c("caseid", "cutoffym")])
origins_national <- fread(here("data/origins_national.csv"))
origins_national <- distinct(origins_national[, c("caseid", "lat_lon_origin_final")])
origins_national <- merge(origins_national, respondents_national)
rm(respondents_national)

# Merge national and local origins data:
origins <- rbindlist(list(origins_local, origins_national), fill = T)
origins[, es_origin := ifelse(is.na(es_origin), 0, es_origin)]
rm(origins_local, origins_national)

# Process origins fields:
origins[, (c("lat", "long")) := list(as.numeric(substr(lat_lon_origin_final, 1, 7)), -as.numeric(substr(lat_lon_origin_final, 10, 16)))]

# Process cutoff date fields:
origins[, (c("year", "mon")) := list(as.numeric(substr(cutoffym, 1, 4)), as.numeric(substr(cutoffym, 5, 6)))]
origins[, (c("lat_lon_origin_final", "cutoffym")) := NULL]

# STEP 3: Create distance matrix ----

# Construct origin spatial object:
origins_st <- st_as_sf(origins, coords = c("long", "lat"),crs="EPSG:4326") %>%
  st_transform(., st_crs("ESRI:102005"))

# Construct destination spatial object:
sites <- sites[order(site_id_agg)]
sites_st <- st_as_sf(sites, coords = c("long", "lat"),crs="EPSG:4326") %>%
  st_transform(., st_crs("ESRI:102005"))

# Calculate O-D distance matrix:
dist <- st_distance(origins_st, sites_st) %>%
  units::set_units(., "mi")
dist <- as.data.table(dist)
names(dist) <- paste0("dist_", 1:ncol(dist))

# Combine distance matrix and origins data:
origins <- cbind(origins, dist)
rm(origins_st, sites, dist)

# STEP 4: Save distance matrix ----

# Save distance matrix:
fwrite(origins, here("data/distance_mat.csv"))

# STEP 5: Figures ----

# Plot distribution of distances:
dist_long <- melt(origins, measure.vars = list(names(origins)[grepl("dist_", names(origins))]), value.name = c("dist"))
respondents_local <- fread(here("data/respondents_local.csv"))
dist_long[, local := ifelse(caseid %in% unique(respondents_local$caseid), "Local Survey", "National Survey")]
fig_dist_freq <- ggplot(dist_long, aes(x = as.numeric(dist), color = local, fill = local, y = after_stat(count))) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme_minimal() +
  theme(plot.subtitle = element_text(hjust = -0.05), legend.title = element_blank(), 
        axis.line.x.bottom = element_line(), axis.line.y.left = element_line(),
        legend.position = c(0.75, 0.75)) +
  xlab(label = expression("One-way Distance (mi.)")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Frequency")
ggsave(here( "output/fig_dist_freq.jpg"), plot = fig_dist_freq, width = 4, height = 4, units = "in", dpi = 500)
rm(origins, dist_long, respondents_local, fig_dist_freq)

# Add site groupings:
site_groups <- data.table(site_id_agg = 1:83)
site_groups[, site_group := ifelse(site_id_agg < 9,
                                   "Texas Sites",
                                   ifelse(
                                     site_id_agg < 35,
                                     "Northern Gulf Sites",
                                     ifelse(site_id_agg < 63, "Peninsula Sites", "South Atlantic Sites")
                                   ))]
sites_st <- merge(sites_st, site_groups, by = "site_id_agg")
rm(site_groups)

# Map sites:
fig_site_map <- ggplot(data = sites_st) +
  annotation_map_tile(type = "osm", zoom = 6, alpha = 0.9) +
  geom_sf(aes(color = site_group, shape = site_group)) +
  theme_bw() +
  scale_color_manual(values = c("#000000", "#0072B2", "#D55E00", "#009E73")) +
  scale_shape_manual(values = c(15, 16, 17, 7)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.275, 0.225), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"))
ggsave(here( "output/fig_site_map.jpg"), plot = fig_site_map, width = 6, height = 4, units = "in", dpi = 500)
rm(fig_site_map, sites_st)
