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

taxa_names <- c('Rhizaria', 'Copepoda',
                'Phaeodaria','Aulacanthidae','Coelodendridae',
                'Aulosphaeridae','Acantharea',
                'Castanellidae','Collodaria',
                'Medusettidae','Foraminifera',
                'Tuscaroridae','living')

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


# remove data.frames with 0's
no_rows <- NULL
for(name in names(rhiz_only$zoo_files)) {
  if(nrow(rhiz_only$zoo_files[[name]]) == 0) {
    no_rows <- c(no_rows, name)
  }
}

rhiz_only$par_files <- rhiz_only$par_files[-which(names(rhiz_only$par_files) %in% no_rows)]
rhiz_only$zoo_files <- rhiz_only$zoo_files[-which(names(rhiz_only$zoo_files) %in% no_rows)]
rhiz_only$meta <- rhiz_only$meta[-which(rhiz_only$meta$profileid %in% no_rows),]

rhiz_sizes <- rhiz_only$zoo_files |> 
  list_to_tib('profileid')
rhiz_sizes$esd <- rhiz_sizes$esd*unique(rhiz_only$meta$acq_pixel)

# |- Density-based filter selection -------------------------
rhiz_sizes$esd |> range()

rhiz_sizes$esd |> 
  density() |> 
  plot()

# |- NBSS option ------------------------------------

calc_nbss <- function(df, metric = 'esd') {
  
  size_limits <- c(min_size = min(df[[metric]]), max_size = max(df[[metric]]))
  
  df$size_class <- df[[metric]] |>
    assign_spectra_bins('log', min_size = size_limits[1], max_size = size_limits[2])
  
  db_num <- df |>
    count_size_classes() |>
    numeric_spectra(vol_L = 1, needs_format = T)
  
  return(db_num)
}


rhiz_nbss <- rhiz_sizes |>
  calc_nbss()

plot(log(n_s,2) ~ mp_size_class, data = rhiz_nbss)



####
# Sample Volume Selection ######################
####

# |- get volume sampled in 25m bins ---------

vol_by_cast <- uvp_data |> 
  get_ecopart_vol() |> 
  lapply(ecopart_vol_bin, depth_breaks = seq(0,1000,25)) |> 
  list_to_tib("profileid") |> 
  bin_format()

vol_by_cast$vol_sampled <- vol_by_cast$vol_sampled * 0.001 # convert to m3

vol_by_cast$zone <- rep(NA, nrow(vol_by_cast))
for(r in 1:nrow(vol_by_cast)) {
  if(vol_by_cast$max_d[r] <=200) {
    vol_by_cast$zone[r] <- 'epi'
  } else if (vol_by_cast$max_d[r] <=1000) {
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

meso_avg <- vol_by_cast$vol_sampled[vol_by_cast$zone == 'meso'] |> 
  mean()

non_detect <- function(prob,vol) {
  log(prob) / -vol
}

non_detect(0.1, epi_avg)
non_detect(0.1, meso_avg)

# |- What about when integrating? ---------------------

# average volumes
epi_avg * (200/25)
meso_avg * ((300)/25)

non_detect(0.1, epi_avg * (200/25))
non_detect(0.1, meso_avg * ((1000-200) /25))


###
# Calculate Binned Abundances #############
###

# |- Rhizaria Only Calcs ----------------------------------------
rhiz_densities <- rhiz_only |> 
  uvp_zoo_conc(breaks = seq(0,1000,25))

# |- Add 0 observations for certain UVP taxa --------------------
all_possible_taxa <- unique(unlist(lapply(rhiz_densities, function(x) x$group)))

# Function to merge each dataframe in the list with the template dataframe

fill_missing_species <- function(df) {
  if(all(all_possible_taxa %in% unique(df$group))) {
    return(df)
  } else {
    merge_df <- expand.grid(
      db = unique(df$db),
      group = all_possible_taxa
    )
    
    df <- df |> 
      merge(merge_df, by = c('group', 'db'), all.y = T)
    
    df[is.na(df)] <- 0
    return(df)
  }
}

corrected_rhiz_densities <- rhiz_densities |>
  lapply(fill_missing_species) |>
  lapply(bin_format)


saveRDS(corrected_rhiz_densities, './data/01_rhiz-densities.rds')
