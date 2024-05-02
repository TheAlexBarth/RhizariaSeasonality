###
# Cast Map
###
library(EcotaxaTools)
library(ggplot2)
library(ggOceanMaps)
library(ggspatial)
library(dplyr)

ctd_data <- readRDS('./data/00_ctd-data.rds')
zoops <- readRDS("./data/00_zoop-uvp.rds")

# need to check ctd vs uvp for correct coords
ctd_coords <- ctd_data |> 
  list_to_tib('cruise_id') |> 
  select(ctd_origfilename, lat, lon, cruise_id, datetime) |> 
  distinct()

coords <- zoops$meta |> 
  select(profileid, latitude, longitude, sampledate,
         ctd_origfilename, tod) |> 
  left_join(ctd_coords, by = "ctd_origfilename")


# |- All cast map ---------------------
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

all_cast_map <- basemap(limits = c(-65,-62,30,33), bathymetry = T, bathy.alpha = 0.5) +  
  geom_point(aes(x = coords$lon, y =coords$lat))+
  theme_bw()+
  labs(x = "",y = "", shape = "")


atlantic_ocean <- ggplot() +
  geom_sf(data = world, fill = '#7b7b7b') +
  geom_rect(aes(xmin = -65, xmax = -62,
                ymin = 30, ymax = 33),
            fill = 'grey', color = 'black', alpha = 0.1) +
  coord_sf(xlim = c(-79,-5), ylim = c(23.5,48))+
  theme_bw()+
  labs(x = "",y = "", shape = "")+
  theme(panel.background = element_rect(fill = 'lightgrey'),
        axis.text.x = element_text(angle = 45, hjust = c(1,1)),
        legend.position = 'bottom')

ggsave('./output/Supplement/s01_background-map.pdf', atlantic_ocean)
ggsave('./output/Supplement/s01_cast-map.pdf', all_cast_map)

