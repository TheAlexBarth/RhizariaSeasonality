###
# Cast Map
###

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

all_cast_map <- ggplot() +
  geom_sf(data = world, fill = 'green')+
  geom_point(aes(x = coords$lon, y =coords$lat))+
  coord_sf(xlim = c(-65,-62), ylim = c(30,33))+
  theme_bw()+
  labs(x = "",y = "", shape = "")+
  theme(panel.background = element_rect(fill = 'cornflowerblue'),
        axis.text.x = element_text(angle = 45, hjust = c(1,1)),
        legend.position = 'bottom')


windows()
