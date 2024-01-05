####
# Import files and format for git repo ###############
####
rm(list = ls())
# # Requires desktop connection

library(EcotaxaTools)
library(batsFtpReadr)
library(suncalc)
library(lubridate)
library(dplyr)
library(ggplot2)
library(sf)

####
# UVP data processing ################
####

# |- Read in -----------------
raw_data <- ecopart_import('~/BATS_data/rhizaria_seasonality/export_raw_20240102_21_52', trim_to_zoo = T)

# |-|- Meta data management -----------------------

# need to edit some stationId by hand
corrected_meta <- read.csv('~/BATS_data/export_all/cast_names_raw.csv') # need to update

# Trim down corrected meta
corrected_meta <- corrected_meta[which(corrected_meta$profileid %in% raw_data$meta$profileid),]
raw_data$meta$stationid <- corrected_meta$stationid
raw_data$meta$ctd_origfilename <- corrected_meta$ctdref
raw_data$meta$programid <- corrected_meta$proj_id




# |-|- Assign Time of Day -----
# apply tod calcs
# There's sometimes I wrote 6 instead of 3 on the latitude
for(i in 1:length(raw_data$meta$latitude)) {
  if(raw_data$meta$latitude[i] > 60) {
    raw_data$meta$latitude[i] <- raw_data$meta$latitude[i] - 30
  }
}

for(i in 1:length(raw_data$meta$longitude)) {
  if(raw_data$meta$longitude[i] > 0) {
    raw_data$meta$longitude[i]  <- -raw_data$meta$longitude[i] 
  }
}

raw_data$meta$tod <- rep(NA, nrow(raw_data$meta))
for(i in 1:length(raw_data$meta$tod)) {
  casttime <- raw_data$meta$sampledate[i]
  suntimes <- getSunlightTimes(
    date = as.Date(casttime),
    lat = raw_data$meta$latitude[i],
    lon = raw_data$meta$longitude[i],
    tz = 'UTC'
  )
  
  if(casttime < suntimes$nauticalDawn | casttime > suntimes$nauticalDusk) {
    raw_data$meta$tod[i] <- 'night'
  } else if (casttime > suntimes$nauticalDawn & casttime < suntimes$nauticalDusk) {
    raw_data$meta$tod[i] <- 'day'
  } else {
    raw_data$meta$tod[i] <- 'twilight'
  }
}

# |- Organizing and removing HS casts -------------------------

# Removing all hs casts then also
hs_index <- which(raw_data$meta$stationid %in% c('HS','hs'))

# need to sort file lists
raw_data$par_files <- raw_data$par_files[order(names(raw_data$par_files), 
                                               raw_data$meta$profileid)]

raw_data$zoo_files <- raw_data$zoo_files[order(names(raw_data$zoo_files),
                                               raw_data$meta$profileid)]

# |-|- slim down the ecopart_obj ------------------------
raw_data$par_files <- raw_data$par_files[-hs_index]
raw_data$zoo_files <- raw_data$zoo_files[-hs_index]
raw_data$meta <- raw_data$meta[-hs_index,]

raw_data <- as_ecopart_obj(raw_data)


# |- Filter out non-living -------------------

zoops <- raw_data |> 
  mod_zoo('names_drop', drop_names = c('not-living', 'darksphere', 'Trichodesmium',
                                       'duplicate', 'othertocheck', 'temp circle', 'house',
                                       'multiple< other', 'part<other', 'aggregate<Radiolaria'),
          drop_children = T)


####
# Import and format CTD data ######################
####

ctd_reader <- function(filename) {
  data <- read.table(filename, header = FALSE, skip = 2)
  
  names(data) <- c('ctd_origfilename', 'dec_year', 'lat', 'lon',
                   'dbar', 'depth', 'temp', 'sal', 'umol_kg',
                   'BAC', 'RFU', 'PAR')
  
  data$cruise_id <- data$ctd_origfilename |> 
    as.character() |> 
    substr(1,5)
  data$lon <- -data$lon
  data$datetime <- date_decimal(data$dec_year)
  
  return(data)
}

ctd_files <- dir('~/BATS_data/rhizaria_seasonality/batsdata/ctd', full.names = T)

ctd_data <- ctd_files |> lapply(ctd_reader)

cruise_ids <- ctd_data |> 
  lapply(`[[`, "cruise_id") |> 
  unlist() |> 
  unique()

names(ctd_data) <- cruise_ids

####
# Geog filters  ################
####

# need to check ctd vs uvp for correct coords
ctd_coords <- ctd_data |> 
  list_to_tib('cruise_id') |> 
  select(ctd_origfilename, lat, lon, cruise_id, datetime) |> 
  distinct()

coords <- zoops$meta |> 
  select(profileid, latitude, longitude, sampledate,
         ctd_origfilename, tod) |> 
  left_join(ctd_coords, by = "ctd_origfilename")
# 
# # |- All cast map ---------------------
# library("rnaturalearth")
# library("rnaturalearthdata")
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# all_cast_map <- ggplot()+
#   geom_sf(data = world, fill = 'green')+
#   geom_point(aes(x = coords$lon, coords$lat,
#                  shape = coords$tod, col = coords$tod))+
#   coord_sf(xlim = c(-65,-62), ylim = c(30,33))+
#   geom_segment(aes(x = -65.1, xend = -62.9,
#                    y = 31.8, yend = 32.75),
#                size = 2, col = 'red')+
#   theme_bw()+
#   labs(x = "",y = "", shape = "")+
#   scale_color_manual(values = c('beige','black'))+
#   theme(panel.background = element_rect(fill = 'cornflowerblue'),
#         axis.text.x = element_text(angle = 45, hjust = c(1,1)),
#         legend.position = 'bottom')

# |- Removing casts outside range ------------------------
drop_casts <- which(
  (coords$lat < 30) |
    (coords$lat > 32.25 & coords$lon < -64)
)

# |- trim and save data ------------------------------------
# note that while creating this data I mess around to get the exact coords
# this code is written to be ran sequentially - be careful to changes

keep_names <- coords$ctd_origfilename[-drop_casts]

trim_zoop <- list()
trim_zoop$meta <- zoops$meta[which(zoops$meta$ctd_origfilename %in% keep_names),]
trim_zoop$par_files <- zoops$par_files[which(names(zoops$par_files) %in% zoops$meta$profileid)]
trim_zoop$zoo_files <- zoops$zoo_files[which(names(zoops$zoo_files) %in% zoops$meta$profileid)]

trim_zoop <- trim_zoop |> as_ecopart_obj()

trim_zoop$meta$cruise_id <- trim_zoop$meta$ctd_origfilename |> 
  as.character() |> 
  substr(1,5) |> 
  unlist() |> 
  as.vector()

saveRDS(trim_zoop, './data/00_zoop-uvp.rds')
saveRDS(ctd_data, './data/00_ctd-data.rds')
saveRDS(coords[-drop_casts,], './data/s00_cast-coords.rds')


####
# Read in bottle data ###################
####

# this needs to be ran after geog filters
bottle_path <- '~/BATS_data/rhizaria_seasonality/batsdata/bats_bottle.txt'
all_bottles <- read.table(bottle_path,
                          header = FALSE, 
                          skip = grep("\\/Variables", readLines(bottle_path)) +1, 
                          fill = TRUE)

# |- Renaming and formatting --------------------
names(all_bottles) <- c('bottle_id',
                        'yyyymmdd',
                        'dec_date',
                        'time',
                        'lat',
                        'lon',
                        'QF',
                        'depth',
                        'temp',
                        'sal_1',
                        'sal_2',
                        'sigth',
                        'O2_umol_kg',
                        'OxFixT',
                        'OxAnom',
                        'CO2',
                        'Alk',
                        'NO3',
                        'NO2',
                        'PO41',
                        'Si',
                        'POC',
                        'PON',
                        'TOC',
                        'TN',
                        'Bact_enumb',
                        'POP',
                        'TDP',
                        "SRP",
                        'BSi',
                        'LSi',
                        'Pro',
                        'Syn',
                        'Piceu',
                        'Naneu')

all_bottles$datetime <- date_decimal(all_bottles$dec_date)
all_bottles$lon <- -all_bottles$lon
all_bottles$ctd_origfilename <- all_bottles$bottle_id |> 
  as.character() |> 
  substr(1,8) |> 
  unlist() |> 
  as.numeric()

# |- Trim and save ----------------
trim_bottles <- all_bottles[which(all_bottles$ctd_origfilename %in% trim_zoop$meta$ctd_origfilename),]

saveRDS(trim_bottles, './data/00_bottles.rds')


####
# Read productivity data #############
####

prod_path <- '~/BATS_data/rhizaria_seasonality/batsdata/bats_primary_production_v003.txt'
prod <- read.table(prod_path,
                          header = FALSE, 
                          skip = grep("\\/Variables", readLines(prod_path)) +1, 
                        fill = TRUE)

# |- rename and format ---------------
names(prod) <- c('bottle_id',
                 'yymmdd_in',
                 'yymmdd_out',
                 'decy_in',
                 'decy_out',
                 'hhmm_in',
                 'hhmm_out',
                 'lat_in',
                 'lat_out',
                 'lon_in',
                 'lon_out',
                 'QF',
                 'depth',
                 'dbar',
                 'temp',
                 'sal',
                 'lt1',
                 'lt2',
                 'lt3',
                 'dark',
                 't0',
                 'pp')


prod$datetime <- date_decimal(prod$decy_in)
prod$lon <- -prod$lon_out
prod$lat <- prod$lat_out
prod$ctd_origfilename <- prod$bottle_id |> 
  as.character() |> 
  substr(1,8) |> 
  unlist() |> 
  as.numeric()

prod$cruise_id <- prod$ctd_origfilename |> 
  as.character() |> 
  substr(1,5) |> 
  unlist() |> 
  as.numeric()

# |- Trim and save ----------------
trim_prod <- prod[which(prod$cruise_id %in% trim_zoop$meta$cruise_id),]

saveRDS(trim_prod, './data/00_prod.rds')


####
# Read in Flux estimates ###############
####

flux_path <- '~/BATS_data/rhizaria_seasonality/batsdata/bats_flux_v003.txt'
flux <- read.table(flux_path,
                   header = FALSE, 
                   skip = grep("\\/Variables", readLines(flux_path)) +1, 
                   fill = TRUE)

names(flux) <- c(
  'cruise_id',
  'depth',
  'yymmdd_in',
  'yymmdd_out',
  'decy_in',
  'decy_out',
  'lat_in',
  'lat_out',
  'lon_in',
  'lon_out',
  'Mass_flux_1',
  'mass_flux_2',
  'mass_flux_3',
  'avg_mass_flux',
  'c_flux_1',
  'c_flux_2',
  'c_flux_3',
  'avg_c_flux',
  'n_flux_1',
  'n_flux_2',
  'n_flux_3',
  'avg_n_flux',
  'p_flux_1',
  'p_flux_2',
  'p_flux_3',
  'avg_p_flux',
  'fbc1','fbc2','fbc3',
  'avg_fbc','fbn1','fbn2',
  'fbn3',
  'avg_fbn'
)

flux$datetime <- date_decimal(flux$decy_in)
flux$lon <- -flux$lon_out
flux$lat <- flux$lat_out


# |- Trim and save ----------------
trim_flux <- flux[which(flux$cruise_id %in% trim_zoop$meta$cruise_id),]

saveRDS(trim_flux, './data/00_flux.rds')
