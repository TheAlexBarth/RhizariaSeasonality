###
# Depth-binning, size-selection, calcs
##


rm(list = ls())
library(EcotaxaTools)

###
# Set up ######################
###

uvp_data <- readRDS("./data/00_zoop-uvp.rds")

####
# Re-naming zoops ########
####

taxa_names <- c('Rhizaria', 'Eumalacostraca', 'Copepoda',
                'Phaeodaria','Aulacanthidae','Coelodendridae',
                'Ostracoda','Aulosphaeridae','Acantharea','Chaetognatha',
                'Castanellidae','Hydrozoa','Collodaria','Actinopterygii',
                'Medusettidae','Foraminifera','Pteropoda',
                'Annelida','Crustacea','living','Decapoda',
                'Tuscaroridae')

uvp_data <- uvp_data |> 
  add_zoo(func = names_to, col_name = 'name',
          new_names = taxa_names)


####
# Size threshold analysis ##########
####

# |- Formatting for Rhiz only ----------------------------------

rhiz_only <- uvp_data |> 
  mod_zoo(func = 'names_keep', keep_names = 'Rhizaria',
          keep_children = T)

# need to account for casts with no observations####

rhiz_sizes <- rhiz_only$zoo_files |> 
  list_to_tib('profileid')



####
# Sample Volume Selection ######################
####

# |- get volume sampled in 10m bins ---------

vol_by_cast <- uvp_data |> 
  get_ecopart_vol() |> 
  lapply(ecopart_vol_bin, depth_breaks = seq(0,1200,20)) |> 
  list_to_tib("profileid") |> 
  bin_format()

vol_by_cast$vol_sampled <- vol_by_cast$vol_sampled * 0.001 # convert to m3

vol_by_cast$zone <- rep(NA, nrow(vol_by_cast))
for(r in 1:nrow(vol_by_cast)) {
  if(vol_by_cast$max_d[r] <=200) {
    vol_by_cast$zone[r] <- 'epi'
  } else if (vol_by_cast$max_d[r] <=1200) {
    vol_by_cast$zone[r] <- 'meso'
  } else {
    vol_by_cast$zone[r] <- 'deep'
  }
}


# |- Detection Probabilities --------------------------
# see barth review and benfield et al 1996 on VPR
# euphotic zone

epi_avg <- vol_by_cast$vol_sampled[vol_by_cast$zone == 'epi'] |> 
  mean()

non_detect <- function(prob,vol) {
  log(prob) / -vol
}

non_detect(0.1, epi_avg)
