rm(list = ls())
graphics.off()


library(ggplot2)
library(dplyr)
library(EcotaxaTools)


# Load the datar

rhiz <- readRDS('./data/00_zoop-uvp.rds') |> 
  mod_zoo(func = names_keep,
          keep_names = 'Rhizaria', keep_children = TRUE)

rhiz$zoo_files <- rhiz$zoo_files[-which(names(rhiz$zoo_files) == 'gf380c11')]
rhiz$par_files <- rhiz$par_files[-which(names(rhiz$par_files) == 'gf380c11')]
rhiz$meta <- rhiz$meta[which(rhiz$meta$profileid != 'gf380c11'),]

rhiz = as_ecopart_obj(rhiz)

vol_s = list()

vol_s <- rhiz |> 
  get_ecopart_vol() |> 
  lapply(
    function(df) {
      sum(df$vol_sampled, na.rm = T)
    }
  ) |> 
  list_to_tib('profileid')



names(vol_s)[1] <- 'volume'


## |- NBSS option ------------------------------------


#' get the nbss for each cast.
calc_nbss <- function(df, metric = 'esd') {
  
  size_limits <- c(min_size = min(df[[metric]]), max_size = max(df[[metric]]))
  
  df$size_class <- df[[metric]] |>
    assign_spectra_bins('log', min_size = size_limits[1], max_size = size_limits[2])
  
  db_num <- df |>
    count_size_classes() |>
    numeric_spectra(vol_L = 1, needs_format = T)
  
  return(db_num)
}


# create bin sizes
rhiz_sizes <- get_all(rhiz, 'esd', pixel_conv = T)
min(rhiz_sizes)
max(rhiz_sizes)


## |-|- Create list of bins




rhiz_nbss <- rhiz_sizes |>
  calc_nbss()

plot(log(n_s,2) ~ mp_size_class, data = rhiz_nbss)
