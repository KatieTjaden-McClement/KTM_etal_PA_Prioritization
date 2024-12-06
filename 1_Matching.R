### PA Priorization
### Tjaden-McClement et al., 2024
### Data Prep and Statistical Matching

list.of.packages <- c("sf",
                      "sp",
                      "rgdal",
                      "rgeos",
                      "tidyverse",
                      "raster",
                      "spdplyr",
                      "spatialEco",
                      "tidyr",
                      "rlist",
                      "purrr",
                      "tibble",
                      "stringr",
                      "pbapply",
                      
                      "ggplot2",
                      "cowplot",
                      "ggsci",
                      "ggmap",
                      "rnaturalearth",
                      
                      "foreach",
                      "doParallel",
                      "parallel",
                      
                      "MatchIt",
                      "Matching",
                      "rgenoud",
                      "CBPS",
                      "optmatch",
                      "cobalt")


lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_classic())

# NOTE: this script uses the "sp", "sf", and "raster" packages, 
# for spatial analysis, which have now been superseded by the 
# "terra" package, which tends to be more efficient for spatial analysis in R

# Also note that certain parts of this script required a 
# virtual machine with 64 GB of ram to run successfully.

##### Prepping World Grid #####

### WDPA Data
# October 2021 WDPA dataset with point data buffered to have their specified area
load("Shapefiles/pa.all.updated100421.Rdata")
pa.all

unique(pa.all$IUCN_CAT)
unique(pa.all$STATUS)

pa.all_sp <- as(pa.all, "Spatial")
pa.all_sp

### Create world grid

## Load in human footprint 1km world grid in Mollweide projection
globe <- raster("Shapefiles/Footprint_v3/wildareas-v3-2009-human-footprint-geotiff/wildareas-v3-2009-human-footprint.tif")
globe

plot(globe)

# Convert grid to data frame - globe is your grid of the world. 
str(globe)

#resample to coarser limit of 5km^2
globe_5 <- raster::aggregate(globe, fact = 5)
plot(globe_5)
globe_5

dat <- as.data.frame(globe_5, xy = TRUE, na.rm = TRUE)

# Convert to spatial object, and convert to same projection
coordinates(dat) <- c("x", "y")

proj4string(dat) <- crs(globe)
str(dat)

# add datID column for merging later
dat$datID <- seq.int(nrow(dat))

#keep just terrestrial grid cells using the Terr Ecoregions shapefile
terr_eco <- st_read("Shapefiles/2. Matching Variables/Terr_Ecoregions/wwf_terr_ecos.shp")
terr_eco <- st_transform(terr_eco, crs(dat))

plot(st_geometry(terr_eco))
str(terr_eco) #has biome, realm, and ecoregion data

terr_eco <- as(terr_eco, "Spatial")

terr_dat <- over(dat, terr_eco)
str(terr_dat)

# merge back with dat to get a spatial points object
terr_dat$datID <- seq.int(nrow(terr_dat))

terr_dat <- sp::merge(dat, terr_dat, 
                      by = "datID",
                      all = T,
                      duplicateGeoms = T)

#need to remove NAs - represent grid cells that aren't terrestrial
terr_dat <- dplyr::filter(terr_dat, ECO_ID != "NA") #5410922 grid cells
plot(terr_dat) #looks good!

dat <- terr_dat

#save for quick access
save(dat, file = "dat_terr_5.RData")

##### Flattening WDPA & Creating Yearly Datasets #####

#remove all columns except the new "ID" column to make the over function run faster
names(pa.all_sp)

pa.all_sp$pa_ID <- seq.int(nrow(pa.all_sp))

pa.all_sp_id <- pa.all_sp[,-(1:26)]
head(pa.all_sp_id)

#dat <- terr_dat_test
str(dat)
dat <- dat %>% 
  dplyr::select(c("datID", "ECO_ID")) #keep only necessary columns
str(dat)

# Frees up memory:
gc(reset = T, full = T)

# needed a VM with 64GB ram to run this line, wouldn't work with 32GB
p_list <- over(dat, pa.all_sp_id, returnList = T) #start 11:45 Oct 20, done by 8pm
str(p_list)

# massive file, save to load in again:
#save(p_list, file = "Key_R_files/p_list.RData")
#load("F:/PA_Prioritization/R files_5/p_list.RData")

# create new dataframe with just year and pa_id columns
glimpse(pa.all)

pa_year <- data_frame("pa_ID" = pa.all_sp$pa_ID, "year" = pa.all_sp$STATUS_YR)
pa_year <- as.data.frame(pa_year)
str(pa_year)

### Converting p_list to yearly datasets
elementID <- list.names(p_list)

# add column to each p_list dataframe with the list element ID = dat grid cell
p_list_ID <- map2(p_list, elementID, ~.x %>% mutate(datID = .y))

# bind all dataframes in p_list together
p_unlisted <- do.call(rbind, p_list_ID)
str(p_unlisted)

# merge back with dat
merge_p_dat <- sp::merge(dat, p_unlisted, 
                         by = "datID",
                         all = T,
                         duplicateGeoms = T)

#merge back with year info
merge_p_year <- sp::merge(merge_p_dat, pa_year,
                          by.x = "pa_ID",
                          all.x = T)

p_allyears <- merge_p_year
head(p_allyears)

# filter to yearly datasets:
for(i in 1980:2021){
  assign(paste0("pa_", i), filter(p_allyears, year <= i))
}

# put yearly datasets in a list for ease of processing
pa_by_year <- mget(ls(pattern = "pa_\\d"))
save(pa_by_year, file = "Key_R_files/pa_by_year.RData")

# flatten yearly datasets, keeping only 1 entry for each datID
pa_by_year_flat <- lapply(pa_by_year, FUN = subset, !duplicated(datID))

#save(pa_by_year_flat, file = "Key_R_files/pa_by_year_flat.RData")
load("F:/PA_prioritization/Key_R_files/pa_by_year_flat.RData")

# create data frame summarizing global PA coverage over time
pa_area_overtime <- lapply(pa_by_year_flat, FUN = length) #why does length work here but not for g200??
pa_area_overtime

# convert to dataframe with year column
pa_over_time <- unlist(pa_area_overtime)
pa_over_time <- as.data.frame(pa_over_time)

str(pa_over_time)

pa_over_time$year <- c(1980:2021)

#go from spatial points represented to actual area covered in km^2 by multiplying
#by 25 km^2 for each point (5x5 grid cell)
pa_over_time$pa_area <- pa_over_time$pa_over_time*25

# plot PA coverage over time
ggplot(pa_over_time, aes(year, pa_area)) +
  geom_point() +
  geom_line()

# plot % PA coverage over time
ggplot(pa_over_time, aes(year, pa_area/(nrow(dat)*25)*100)) +
  geom_point(size = 4) +
  geom_line(size = 1.5) +
  labs(x = "Year", y = "% terrestrial area protected") +
  ylim(0,20) +
  theme(text = element_text(size = 20))

##### Biodiversity hotspots #####
hotspots <- st_read("Shapefiles/1. Prioritization_schemes/Biodiversity Hotspots_old/data/hotspots_revisited_2004_polygons.shp")
hotspots

plot(st_geometry(hotspots)) #looks good
crs(hotspots) #in WGS84

hotspots_terr <- filter(hotspots, TYPE == "hotspot_area") #removes outer limit polygons (marine)

hotspots_sp <- as(hotspots_terr, "Spatial")

# reprojection in R was causing weird distortions... did in ArcGIS
# but might work better using newer terra package

# write to reproject in ArcGIS:
# writeOGR(hotspots_sp,layer = "hotspots_sp",
#         "hotspots_sp", driver = "ESRI Shapefile")

# reprojected in ArcGIS
hotspots_moll <- st_read("Shapefiles/1. Prioritization_schemes/Hotspots/hotspots_moll/hotspots_moll.shp")
crs(hotspots_moll)

hotspots_map <- ggplot(countries_moll) +
  geom_sf(fill = "#EFEFEF",
          color = "#D0D0D0",
          size = 0.1) +
  geom_sf(data = hotspots_moll,
          aes(colour = "#1F77B4FF"),
          fill = "#1F77B4FF") +
  scale_colour_identity(name = "",
                        labels = "Hotspots",
                        breaks = "#1F77B4FF",
                        guide = "legend") +
  theme(legend.position = c(0.9, 0.9),
        legend.background =  element_rect(colour = NA))

hotspots_sp <- as(hotspots_moll, "Spatial")

plot(hotspots_sp) # looks good

hotspots_over <- over(dat, hotspots_sp)
#save(hotspots_over, file = "hotspots.RData")

### Protection in Hotspots

#load("F:/PAs_R/hotspots_5.RData")
hotspots <- hotspots_over

# add datID to link to dat
hotspots$datID <- elementID
str(hotspots)

hotspots <- sp::merge(dat, hotspots, 
                      by = "datID",
                      all = T,
                      duplicateGeoms = T)

#need to separate out hotspots vs. not hotspot grid cells
not_hotspots <- filter(hotspots, is.na(NAME))
#plot(not_hotspots) #looks good
#save(not_hotspots, file = "not_hotspots_5.RData")

hotspots <- filter(hotspots, NAME != "NA")
#plot(hotspots)
#save(hotspots, file = "hotspots_5.RData")

## Find overlap between hotspots and protected areas
# NOTE: this does not track hotspot cells that aren't protected as needed for panel data
protected_hotspots <- lapply(pa_by_year_flat, FUN = filter, datID %in% hotspots$datID)

pa_hotspots_overtime <- lapply(protected_hotspots, FUN = nrow) 
pa_hotspots_overtime <- unlist(pa_hotspots_overtime) 
pa_hotspots_overtime <- as.data.frame(pa_hotspots_overtime)

# add column to pa_over_time
pa_over_time <- cbind(pa_over_time, pa_hotspots_overtime)

# convert to area
pa_over_time$hotspots_coverage_area <- pa_over_time$pa_hotspots_overtime*25

# Plotting protected area growth inside vs outside of hotspots from 1980-2019:
ggplot(data = pa_over_time) +
  geom_point(aes(x = year, y = hotspots_coverage_area, 
                 colour = "Biodiversity Hotspots"), 
             size = 2.5, show.legend = T) +
  geom_line(aes(x = year, y = hotspots_coverage_area),
            color = "#1F77B4FF", size = 1) +
  geom_point(aes(x = year, y = pa_area - hotspots_coverage_area,
                 colour = "Outside Biodiversity hotspots"),
             size = 2.5) +
  geom_line(aes(x = year, y = pa_area - hotspots_coverage_area),
            size = 1) +
  scale_color_manual(values = c("#1F77B4FF", "black")) +
  labs(y = "Protected area (km^2)", x = "Year", colour = "")

##### Last of the Wild #####

lotw_raw <- st_read("Shapefiles/1. Prioritization_schemes/LOTW/ltw-global-geo/ltw_v2geo.shp")
lotw_raw
crs(lotw_raw) #CRS arguments: +proj=longlat +ellps=GRS80 +no_defs 

plot(st_geometry(lotw_raw))

# this dataset contains all the 10% most wild areas of each biome of each realm,
# not the 10 largest polygons of the 10% wildest area
# need to keep only the 10 largest polygons (excluding any less than 5km^2) in each biome and realm to get to LOTW

lotw_raw <- st_transform(lotw_raw, crs(dat))
lotw_raw$area <- st_area(lotw_raw)
lotw_raw$area <- as.numeric(lotw_raw$area)

# Confirm reporjection didn't make anything go crazy
plot(st_geometry(lotw_raw)) #looks good

lotw <- lotw_raw %>% 
  filter(!is.na(REALM), area > 5e+6) %>% 
  group_by(REALM, BIOME) %>% 
  arrange(desc(area), .by_group = T) %>% 
  slice(1:10)
# results in 567 sites

plot(st_geometry(lotw))

lotw_map <- ggplot(countries_moll) +
  geom_sf(fill = "#EFEFEF",
          color = "#D0D0D0",
          size = 0.1) +
  geom_sf(data = lotw,
          fill = "#9769C1",
          aes(colour = "#9769C1")) +
  scale_colour_identity(name = "",
                        labels = "LOTW",
                        breaks = "#9769C1",
                        guide = "legend") +
  theme(legend.position = c(0.9, 0.9),
        legend.background =  element_rect(colour = NA))

# now looks same as LW figure in Brooks et al (2006)

lotw_sp <- as(lotw, "Spatial")
lotw_over <- over(dat, lotw_sp)

lotw_over$datID <- dat$datID #dat is terrestrial subset, datID doesn't start at 1

# seperate out lotw cells from non-lotw cells
lotw_datIDs <- filter(lotw_over, !is.na(ID)) %>% 
  pull(datID) # 1397574

not_lotw_datIDs <- filter(lotw_over, is.na(ID)) %>% 
  pull(datID) # 4013348

lotw <- dat %>% 
  filter(datID %in% lotw_datIDs)
#plot(lotw)

not_lotw <- dat %>% 
  filter(datID %in% not_lotw_datIDs)
#plot(not_lotw)

protected_lotw <- lapply(pa_by_year_flat, FUN = filter, datID %in% lotw$datID)

pa_lotw_overtime <- lapply(protected_lotw, FUN = nrow)
pa_lotw_overtime <- unlist(pa_lotw_overtime) 
pa_lotw_overtime <- as.data.frame(pa_lotw_overtime)

pa_lotw_overtime$year <- c(1980:2021)

#repeat for cells that aren't in lotw
protected_not_lotw <- lapply(pa_by_year_flat, FUN = filter, datID %in% not_lotw$datID)

pa_not_lotw_overtime <- lapply(protected_not_lotw, FUN = nrow)
pa_not_lotw_overtime <- unlist(pa_not_lotw_overtime)

pa_lotw_overtime$pa_not_lotw_overtime <- pa_not_lotw_overtime

pa_lotw_overtime <- pa_lotw_overtime %>% 
  mutate(pa_lotw_area = pa_lotw_overtime*25,
         pa_not_lotw_area = pa_not_lotw_overtime*25,
         prop_lotw = pa_lotw_area/(nrow(lotw)*25),
         prop_not_lotw = pa_not_lotw_area/(nrow(not_lotw)*25))
head(pa_lotw_overtime)

#Plot area protected inside vs. outside LOTW sites over time
ggplot(pa_lotw_overtime) +
  geom_point(aes(x = year, y = pa_lotw_area,
                 colour = "Last of the Wild"), 
             size = 3.5, show.legend = T) +
  geom_line(aes(x = year, y = pa_lotw_area,
                colour = "Last of the Wild"),
            size = 1.5) +
  geom_point(aes(x = year, y = pa_not_lotw_area,
                 colour = "Outside Last of the Wild"),
             size = 3.5) +
  geom_line(aes(x = year, y = pa_not_lotw_area),
            size = 1.5) +
  scale_color_manual(values = c("#9769C1", "black")) +
  labs(y = "Area protected (km^2)", x = "Year", color = "")

# Plot PROPORTION area protected inside vs. outside LOTW
ggplot(pa_lotw_overtime) +
  geom_point(aes(x = year, y = prop_lotw,
                 colour = "Last of the Wild"), 
             size = 3.5, show.legend = T) +
  geom_line(aes(x = year, y = prop_lotw,
                colour = "Last of the Wild"),
            size = 1.5) +
  geom_point(aes(x = year, y = prop_not_lotw,
                 colour = "Outside Last of the Wild"),
             size = 3.5) +
  geom_line(aes(x = year, y = prop_not_lotw),
            size = 1.5) +
  scale_color_manual(values = c("#9769C1", "black")) +
  labs(y = "Proportion area protected", x = "Year", color = "")

##### Matching Variables #####
#frees up memory:
gc(reset = T, full = T)

### Load Matching Variables

# Human footprint v3 in Moll aggregated to 5km^2 cells as base to reproject others
# (repeat of early step)
globe <- raster("Shapefiles/Human Footprint/wildareas-v3-2009-human-footprint.tif")
globe_5 <- raster::aggregate(globe, fact = 5)

## Human Footprint v2
footprint_raw <- raster("Shapefiles/2. Matching Variables/Human Footprint/hfp_global_geo_grid/hf_v2geo")
plot(footprint_raw) # in ellps=clrk66 projection

footprint <- projectRaster(from = footprint_raw, to = globe_5, method = "bilinear")
footprint
plot(footprint)

## Elevation
elevation_raw <- raster("Shapefiles/2. Matching Variables/Elevation/mn30_grd")
elevation_raw #in WGS84, resolution: 0.008333333, 0.008333333 (x, y)
plot(elevation_raw)

# reproject to Mollwiede
elevation <- projectRaster(from = elevation_raw, to = globe_5, method = "bilinear") #gave a bunch of warnings but I think it's okay
elevation
plot(elevation) #looks good

## Agricultural potential
ag_opp_raw <- raster("Shapefiles/2. Matching Variables/Ag. Potential/argmean_max/agrmean_max")
plot(ag_opp_raw)

ag_opp <- projectRaster(from = ag_opp_raw, to = globe_5, method = "bilinear")
ag_opp #worked - same projection and resolution
plot(ag_opp)

#save(ag_opp_moll, file = "ag_opp_moll.RData")

## Population Density
pop_den_raw <- raster("Shapefiles/2. Matching Variables/Population Density/gluds00ag")#adjusted
# *ag files are adjusted to match UN country totals
pop_den_raw
plot(pop_den_raw)

pop_den <- projectRaster(from = pop_den_raw, to = globe_5, method = "bilinear")
pop_den
plot(pop_den)

save(pop_den, file = "Key_R_files/Matching/pop_den.RData")

## Global roads
roads_raw <- raster("Shapefiles/2. Matching Variables/Roads_new/Roads.tif")
roads_raw
plot(roads_raw) #already in Mollwiede, but has 1km^2 res

roads <- projectRaster(from = roads_raw, to = globe_5, method = "bilinear")
plot(roads)

### Extract values to dat

#save(elevation, file = "elevation.RData")
#save(footprint, file  = "footrpint.RData")

matching_vars <- stack(elevation, footprint, ag_opp, pop_den, roads)

#Extract values of matching vars across dat spatial points
matching_values <- raster::extract(matching_vars, dat)

summary(matching_values) #lots of NAs, especially in ag_opp: 208565 - only 4% of terr dat cells though, ok

dat_matching <- cbind(dat, matching_values)
str(dat_matching)

##### Matching Hotspots #####
#load("Key_R_files/hotspots_5.RData")
str(hotspots)
hotspots$hotspot <- 1

#load("Key_R_files/not_hotspots_5.RData")
str(not_hotspots)
not_hotspots$hotspot <- 0

dat_hotspots <- rbind(hotspots, not_hotspots)
head(dat_hotspots)

# Keep just datID and whether it's a hotspot
datID_hotspot <- dplyr::select(dat_hotspots, "datID", "hotspot")
head(datID_hotspot)

# merge with matching covarite data for dat cells
dat_matching_hotspots <- merge(dat_matching, datID_hotspot, by = "datID")
head(dat_matching_hotspots) #this has matching variable info and whether a dat cell is a hotspot!

#Can't have any missing values in the covariates:
#matching_data <- sp.na.omit(dat_matching_hotspots) # causing R to abort
summary(matching_data)
str(matching_data)

matching_data <- dat_matching_hotspots %>% 
  dplyr::select(datID, REALM, BIOME, ECO_ID, eco_code, mn30_grd,
                hf_v2geo, agrmean_max, gluds00ag, Roads, hotspot, x, y) %>% 
  rename(elevation = mn30_grd,
         footprint = hf_v2geo,
         ag = agrmean_max,
         pop = gluds00ag,
         roads = Roads) %>% 
  filter(!is.na(elevation),
         !is.na(footprint),
         !is.na(ag),
         !is.na(pop),
         !is.na(roads)) #removes 240,666 cells

summary(matching_data)

# add country info to dat
countries_raw <- st_read("Shapefiles/2. Matching Variables/Countries/countries.shp")

#reproject to mollweide
countries_moll <- st_transform(countries_raw, crs(dat))
plot(st_geometry(countries_moll))

countries <- as(countries_moll, "Spatial")

dat_countries <- over(dat, countries)

# filter to just keep country and datID 
head(dat_countries)

datID_country <- dat %>% 
  mutate(country = dat_countries$COUNTRY) %>% 
  dplyr::select(datID, country)

datID_country <- as.data.frame(datID_country) %>% 
  dplyr::select(-x, -y)
head(datID_country)

# merge with matching data by datID
matching_data <- merge(matching_data, datID_country, by = "datID", all.x = T) %>% 
  filter(!is.na(country))
head(matching_data)

matching_data_grouped_country <- matching_data %>% 
  as.data.frame() %>% 
  group_by(country) %>% 
  group_split() %>% 
  as.list()
# only 168 country groups

length(unique(countries_moll$COUNTRY)) #249 countries
length(unique(dat_countries$COUNTRY)) #235 countries
length(unique(matching_data$country)) #168 countries

setdiff(countries_moll$COUNTRY, dat_countries$COUNTRY) #all tiny countries e.g. vatican city
setdiff(countries_moll$COUNTRY, matching_data$country) #almost all islands and tiny countries...
# We are losing these countries because one of the other layers is missing data for their grid cells

### CBPS matching across country & biome groups - with calipers
head(matching_data)

matching_groups_cb <- matching_data %>% 
  as.data.frame() %>% 
  group_by(country, BIOME) %>% 
  group_split() %>% 
  as.list()
# 584 groups

# Keep only country-biome groups with both hotspot and non-hotspot grid cells
matching_groups_cb <- list.filter(matching_groups_cb, length(unique(hotspot)) == 2)
# 211 groups

# do matching as loop first to get rid of any groups that don't have any matches because of the calipers:
# for(i in 208:length(matching_groups_cb)){
#   assign(paste0("hotspot_cb_match_group_test_", i),
#          matchit(formula = hotspot ~ elevation +
#                    footprint + ag + pop + roads,
#                  data = matching_groups_cb[[i]],
#                  method = "nearest",
#                  distance = "cbps",
#                  caliper = c(elevation = 0.25, footprint = 0.25,
#                              ag = 0.25, pop = 0.25, roads = 0.25),
#                  std.caliper = c(TRUE, TRUE, TRUE, TRUE, TRUE)))
# 
#   print(paste(i, "group(s) matched"))
# }
# Errors on groups 1, 2, 5, 6, 17, 18, 48, 52, 59, 69, 81, 91, 94, 99, 104, 109, 110, 113, 116, 125
# 134, 142, 151, 160, 169, 180, 181, 183, 188, 190, 205, 207

# filter out those groups
matching_groups_cb_cal <- matching_groups_cb[-c(1, 2, 5, 6, 17, 18, 48, 52, 59, 69, 81, 91, 94, 
                                                99, 104, 109, 110, 113, 116, 125, 134, 142, 151, 
                                                160, 169, 180, 181, 183, 188, 190, 205, 207)]
# left with 179 groups

# run the matching across those groups
hotspot_cb_cal_match_list <- pblapply(X = matching_groups_cb_cal, 
                                      FUN = matchit,
                                      formula = hotspot ~ elevation +
                                        footprint + ag + pop + roads,
                                      method = "nearest",
                                      distance = "cbps",
                                      caliper = c(elevation = 0.25, footprint = 0.25,
                                                  ag = 0.25, pop = 0.25, roads = 0.25),
                                      std.caliper = c(TRUE, TRUE, TRUE, TRUE, TRUE))

for(i in 1:length(hotspot_cb_cal_match_list)){
  assign(paste0("hotspot_cb_cal_match_group_", i), 
         match.data(object = hotspot_cb_cal_match_list[[i]], 
                    data = matching_groups_cb_cal[[i]]))
}


list_hotspot_cb_cal <- mget(ls(pattern = "hotspot_cb_cal_match_group_"))
list_hotspot_cb_cal <- lapply(list_hotspot_cb_cal, as_tibble)

hotspot_match_data_cb_cal <- do.call(rbind, list_hotspot_cb_cal)
save(hotspot_match_data_cb_cal, file = "Key_R_files/Matching/Hotspots/hotspot_match_data_cb_cal.RData")

# look at balance after matching
hotspot_cb_cal_match <- matchit(formula = hotspot ~ elevation +
                                  footprint + ag + pop + roads, 
                                data = hotspot_match_data_cb_cal,
                                method = NULL,
                                distance = "glm")
summary(hotspot_cb_cal_match) # n = 128,392

var_names <- data.frame(old = c("elevation", "ag", "footprint",
                                "pop", "roads"),
                        new = c("Elevation", "Ag. potential",
                                "Human footprint", "Pop. dens",
                                "Road dens."))

# Love plot to look at matching effectiveness
hotspot_love <- love.plot(hotspot_cb_cal_match, colours = "#1F77B4FF",
                          var.names = var_names,
                          sample.names = "Matched Hotspots") # balance very good!!

# Make table of protected area over time in control and treated cells
hotspot_cb_cal_overtime <- DiD_table(matched_data = hotspot_match_data_cb_cal, 
                                     ref_data = pa_by_year_flat,
                                     ID_col = datID, 
                                     trt = hotspot,
                                     years = c(1980:2021),
                                     grid_cell_area = 25)

# plot trends in protection over time
DiD_plot_prop(data = hotspot_cb_cal_overtime,
              x = year,
              trt_column = trt_prop,
              con_column = cons_prop,
              year_trt_est = 2000,
              trt_name = "Biodiversity Hotpsots",
              trt_colour = "#1F77B4FF")

cb_cal_match_mapping <- dat %>% 
  filter(datID %in% c(hotspot_match_data_cb_cal$datID)) %>% 
  mutate(treatment = ifelse(datID %in% hotspots$datID,
                            "hotspot",
                            "control"))

cb_cal_match_mapping <- st_as_sf(cb_cal_match_mapping)

hotspot_cb_cal_match_map <- ggplot(countries_moll) +
  geom_sf(fill = "#EFEFEF",
          color = "#D0D0D0") +
  geom_sf(data = cb_cal_match_mapping,
          aes(color = treatment),
          size = 0.1,
          alpha = 0.05) +
  scale_fill_manual(values = c("#595757", "#1F77B4FF"), 
                    aesthetics = "color",
                    guide = guide_legend(override.aes =
                                           list(fill = c("#595757",
                                                         "#1F77B4FF"),
                                                size = 3,
                                                alpha = 1)),
                    labels = c("Matched Controls",
                               "Matched Hotspots")) +
  labs(color = "") +
  theme(legend.position = c(0.9, 0.75),
        legend.background = element_rect(fill = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

## Visualize with only grid cells that weren't protected in 1980
pa_1980 <- as.data.frame(pa_by_year_flat["pa_1980"])

hotspot_cb_cal_unprotected <- hotspot_match_data_cb_cal %>% 
  filter(!(datID %in% pa_1980$pa_1980.datID))

# Make table of protected area over time in control and treated cells
hotspot_cb_cal_unprotected_overtime <- DiD_table(matched_data = hotspot_cb_cal_unprotected, 
                                                 ref_data = pa_by_year_flat,
                                                 ID_col = datID, 
                                                 trt = hotspot,
                                                 years = c(1980:2021),
                                                 grid_cell_area = 25)
str(hotspot_cb_cal_unprotected_overtime)

# plot trends in protection over time
DiD_plot(data = hotspot_cb_cal_unprotected_overtime,
         time_var = year,
         trt_col = trt_area,
         con_col = cons_area,
         time_trt_est = 2000,
         trt_name = "Biodiversity Hotpsots",
         trt_colour = "#1F77B4FF",
         ylab = "Area protected (km^2)")

##### Matching LOTW #####

# add row with whether the cell is in LOTW or not
lotw_datIDs <- as.vector(lotw$datID)

matching_data <- matching_data %>% 
  dplyr::mutate(lotw = if_else(datID %in% lotw_datIDs,
                               1, 0))
head(matching_data)

# Check pre-match imbalance
pre_match_lotw <- matchit(lotw ~ elevation +
                            footprint + ag + pop + roads,
                          data = matching_data,
                          method = NULL,
                          distance = "glm")
summary(pre_match_lotw)
love.plot(pre_match_lotw)
# MASSIVE difference in population density

lotw_matching_groups_cb <- list.filter(matching_groups_cb, length(unique(lotw)) == 2)
# 186 groups

# do as loop first to get rid of any groups that don't have any matches because of the calipers:
# for(i in 164:length(lotw_matching_groups_cb)){
#   assign(paste0("lotw_cb_match_group_test_", i),
#          matchit(formula = lotw ~ elevation +
#                    footprint + ag + pop + roads,
#                  data = lotw_matching_groups_cb[[i]],
#                  method = "nearest",
#                  distance = "cbps",
#                  caliper = c(elevation = 0.25, footprint = 0.25,
#                              ag = 0.25, pop = 0.25, roads = 0.25),
#                  std.caliper = c(TRUE, TRUE, TRUE, TRUE, TRUE)))
# 
#   print(paste(i, "group(s) matched"))
# }
# Errors on groups 10, 17, 18, 58, 74, 101, 123, 134, 138, 140, 146, 163

# [1] "73 group(s) matched"
# Error: Calipers cannot be used with binary, factor, or character variables. Offending variables:
#   roads
group74 <- as.data.frame(lotw_matching_groups_cb[74])
# value for roads is 0 for everything, remove this group too

lotw_matching_groups_cb_cal <- lotw_matching_groups_cb[-c(10, 17, 18, 58, 74, 101, 123, 134, 138, 140, 146, 163)]

lotw_match_cb_cal <- pblapply(X = lotw_matching_groups_cb_cal, 
                              FUN = matchit,
                              formula = lotw ~ elevation +
                                footprint + ag + pop + roads,
                              method = "nearest",
                              distance = "cbps",
                              caliper = c(elevation = 0.25, footprint = 0.25,
                                          ag = 0.25, pop = 0.25, roads = 0.25),
                              std.caliper = c(TRUE, TRUE, TRUE, TRUE, TRUE))
# took 5 hours

# get match data
for(i in 1:length(lotw_match_cb_cal)){
  assign(paste0("lotw_cb_cal_match_group_", i), 
         match.data(object = lotw_match_cb_cal[[i]], 
                    data = lotw_matching_groups_cb_cal[[i]]))
}


list_cb_cal <- mget(ls(pattern = "lotw_cb_cal_match_group_"))
list_cb_cal <- lapply(list_cb_cal, as_tibble)

lotw_match_data_cb_cal <- do.call(rbind, list_cb_cal)
#save(lotw_match_data_cb_cal, file = "Key_R_files/Matching/LOTW/lotw_match_data_cb_cal.RData")


#look at overall balance after matching
lotw_cb_cal_match <- matchit(formula = lotw ~ elevation +
                               footprint + ag + pop + roads, 
                             data = lotw_match_data_cb_cal,
                             method = NULL,
                             distance = "glm")
summary(lotw_cb_cal_match) # n = 484,489

var_names <- data.frame(old = c("elevation", "ag", "footprint",
                                "pop", "roads"),
                        new = c("Elevation", "Ag. potential",
                                "Human footprint", "Pop. dens",
                                "Road dens."))

lotw_love <- love.plot(lotw_cb_cal_match, colours = "#9769C1",
                       var.names = var_names,
                       sample.names = "Matched LOTW",
                       title = "") # pretty good match, not as good as hotspots

lotw_cb_cal_overtime <- DiD_table(matched_data = lotw_match_data_cb_cal, 
                                  ref_data = pa_by_year_flat,
                                  ID_col = datID, 
                                  trt = lotw,
                                  years = c(1980:2021),
                                  grid_cell_area = 25)

DiD_plot(data = lotw_cb_cal_overtime,
         time_var = year,
         trt_col = trt_prop,
         con_col = cons_prop,
         time_trt_est = 2002,
         trt_name = "Last of the Wild",
         trt_colour = "#9769C1",
         ylab = "Proportion area protected")

lotw_cb_cal_match_mapping <- dat %>% 
  filter(datID %in% c(lotw_match_data_cb_cal$datID)) %>% 
  mutate(treatment = ifelse(datID %in% lotw$datID,
                            "lotw",
                            "control"))

lotw_cb_cal_match_mapping <- st_as_sf(lotw_cb_cal_match_mapping)

lotw_cb_cal_match_map <- ggplot(countries_moll) +
  geom_sf(fill = "#EFEFEF",
          color = "#D0D0D0") +
  geom_sf(data = lotw_cb_cal_match_mapping,
          aes(color = treatment),
          size = 0.1,
          alpha = 0.05) +
  scale_fill_manual(values = c("#595757", "#9769C1"), 
                    aesthetics = "color",
                    guide = guide_legend(override.aes =
                                           list(fill = c("#595757",
                                                         "#9769C1"),
                                                size = 3,
                                                alpha = 1)),
                    labels = c("Matched Controls",
                               "Matched LOTW")) +
  labs(color = "") +
  theme(legend.position = c(0.9, 0.75),
        legend.background = element_rect(fill = NA))

## Visualize with only grid cells that weren't protected in 1980
pa_1980 <- as.data.frame(pa_by_year_flat["pa_1980"])

lotw_cb_cal_unprotected <- lotw_match_data_cb_cal %>% 
  filter(!(datID %in% pa_1980$pa_1980.datID))

lotw_cb_cal_unprotected_overtime <- DiD_table(matched_data = lotw_cb_cal_unprotected, 
                                              ref_data = pa_by_year_flat,
                                              ID_col = datID, 
                                              trt = lotw,
                                              years = c(1980:2021),
                                              grid_cell_area = 25)

DiD_plot(data = lotw_cb_cal_unprotected_overtime,
         time_var = year,
         trt_col = trt_prop,
         con_col = cons_prop,
         time_trt_est = 2002,
         trt_name = "Last of the Wild",
         trt_colour = "#9769C1",
         ylab = "Proportion area protected")

### Plot distribution of matching covariates in total sample vs matched

cov_plotting <- matching_data %>% 
  as_tibble() %>%
  mutate(lotw = as.factor(lotw),
         hotspot = as.factor(hotspot),
         log_ag = log(ag+1),
         log_pop = log(pop+1),
         log_roads = log(roads+1),
         log_footprint = log(footprint+1)) %>% 
  dplyr::select(datID, elevation, log_footprint, log_ag, log_pop, log_roads, 
                lotw, hotspot) %>% 
  pivot_longer(cols = c(elevation, log_footprint, log_ag, log_pop, log_roads),
               names_to = "covariate",
               values_to = "cov_value") %>% 
  mutate(covariate = factor(covariate, 
                            labels = c("Elevation (m)", "log(Ag. potential)",
                                       "log(Human footprint)", "log(Pop. dens. (people/km^2))",
                                       "log(Road dens. (km/km^2))")))

# covariate values inside vs. outside LOTW
LOTW_pre_match_covs <- ggplot(cov_plotting, 
                              aes(x = cov_value, group = lotw, 
                                  fill = lotw)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#595757", "#9769C1"),
                    aesthetics = "fill",
                    labels = c("Outside LOTW", "LOTW")) +
  facet_wrap(~covariate, scales = "free",
             labeller = label_value) +
  labs(y = "", x = "", fill = "") +
  theme(legend.position = c(0.85, 0.25))

# covariate values inside vs. outside Hotspots
hotspots_pre_match_covs <- ggplot(cov_plotting, 
                                  aes(x = cov_value, group = hotspot, 
                                      fill = hotspot)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#595757", "#1F77B4FF"),
                    aesthetics = "fill",
                    labels = c("Outside Hotspots", "Hotspots")) +
  facet_wrap(~covariate, scales = "free",
             labeller = label_value) +
  labs(y = "", x = "", fill = "") +
  theme(legend.position = c(0.85, 0.25))

## Covariate balance after matching for LOTW
lotw_match_cov_plot <- matching_data %>% 
  as_tibble() %>% 
  filter(datID %in% lotw_match_data_cb_cal$datID) %>% 
  mutate(lotw = as.factor(lotw),
         log_ag = log(ag+1),
         log_pop = log(pop+1),
         log_roads = log(roads+1),
         log_footprint = log(footprint+1)) %>% 
  dplyr::select(datID, elevation, log_footprint, log_ag, log_pop, log_roads, 
                lotw) %>% 
  pivot_longer(cols = c(elevation, log_footprint, log_ag, log_pop, log_roads),
               names_to = "covariate",
               values_to = "cov_value") %>% 
  mutate(covariate = factor(covariate, 
                            labels = c("Elevation (m)", "log(Ag. potential)",
                                       "log(Human footprint)", "log(Pop. dens. (people/km^2))",
                                       "log(Road dens. (km/km^2))")))

lotw_matched_covs <- ggplot(lotw_match_cov_plot, 
                            aes(x = cov_value, group = lotw, 
                                fill = lotw)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#595757", "#9769C1"),
                    aesthetics = "fill",
                    labels = c("Matched Controls", "Matched LOTW")) +
  facet_wrap(~covariate, scales = "free",
             labeller = label_value) +
  labs(y = "", x = "", fill = "") +
  theme(legend.position = c(0.85, 0.25))

## Covariate balance after matching for Hotspots
hotspot_match_cov_plot <- matching_data %>% 
  as_tibble() %>% 
  filter(datID %in% hotspot_match_data_cb_cal$datID) %>% 
  mutate(hotspot = as.factor(hotspot),
         log_ag = log(ag+1),
         log_pop = log(pop+1),
         log_roads = log(roads+1),
         log_footprint = log(footprint+1)) %>% 
  dplyr::select(datID, elevation, log_footprint, log_ag, log_pop, log_roads, 
                hotspot) %>% 
  pivot_longer(cols = c(elevation, log_footprint, log_ag, log_pop, log_roads),
               names_to = "covariate",
               values_to = "cov_value") %>% 
  mutate(covariate = factor(covariate, 
                            labels = c("Elevation (m)", "log(Ag. potential)",
                                       "log(Human footprint)", "log(Pop. dens. (people/km^2))",
                                       "log(Road dens. (km/km^2))")))

hotspots_matched_covs <- ggplot(hotspot_match_cov_plot, 
                                aes(x = cov_value, group = hotspot, 
                                    fill = hotspot)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#595757", "#1F77B4FF"),
                    aesthetics = "fill",
                    labels = c("Matched Controls", "Matched Hotspots")) +
  facet_wrap(~covariate, scales = "free",
             labeller = label_value) +
  labs(y = "", x = "", fill = "") +
  theme(legend.position = c(0.85, 0.25))

matching_covs_before_after <- cowplot::plot_grid(hotspots_pre_match_covs, LOTW_pre_match_covs,
                                                 hotspots_matched_covs, lotw_matched_covs,
                                                 labels = "auto")
ggsave("matching_covs_before_after.png", width = 13, height = 8)

# love plots together
cb_cal_love_plots <- plot_grid(hotspot_love, #from Hotspot_DiD script
                               lotw_love,
                               labels = "auto",
                               nrow = 2)
ggsave("cb_cal_love_plots.png", width = 8, height = 6)

# nicer maps of schemes together:
scheme_maps <- plot_grid(hotspots_map, lotw_map,
                         labels = "auto",
                         nrow = 2)
ggsave("scheme_maps.png", width = 6, height = 8)

# nicer maps of matched cb_cal samples together
matched_maps <- plot_grid(hotspot_cb_cal_match_map,
                          lotw_cb_cal_match_map,
                          labels = "auto",
                          nrow = 2)
ggsave("matched_sample_maps.png", width = 9, height = 12)
