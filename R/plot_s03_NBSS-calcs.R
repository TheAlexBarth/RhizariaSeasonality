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
calc_nbss_width <- function(df,
                      metric = 'esd',
                      conv = 0.092,
                      vol) {
  
  # assign size bins
  # using width with prefixed sizes will make all bins the same regardless
  # of the df input
  size_classes <- (df[[metric]] * conv) |>
    assign_spectra_bins('width',
                        width = 1,
                        min_size = size_limits[1], 
                        max_size = size_limits[2])
  
  ns = (table(size_classes) / vol)
  
  ns[which(ns == 0)] <- NA
  
  outdf <- data.frame(
    size_bin = levels(size_classes),
    count = as.vector(table(size_classes)),
    ns = as.vector(ns)
  ) |> 
    bin_format(col = 'size_bin')
  return(outdf)
}


calc_nbss_log <- function(df,
                            metric = 'esd',
                            conv = 0.092,
                            vol) {
  
  # assign size bins
  # using width with prefixed sizes will make all bins the same regardless
  # of the df input
  size_classes <- (df[[metric]] * conv) |>
    assign_spectra_bins('log',
                        min_size = size_limits[1], 
                        max_size = size_limits[2])
  
  ns = (table(size_classes) / vol) |> 
    log()
  
  ns[is.infinite(ns)] <- NA
  
  outdf <- data.frame(
    size_bin = levels(size_classes),
    count = as.vector(table(size_classes)),
    ns = as.vector(ns)
  ) |> 
    bin_format(col = 'size_bin')
  return(outdf)
}



# create bin sizes
rhiz_sizes <- get_all(rhiz, 'esd', pixel_conv = T)
size_limits = c(min(rhiz_sizes), max(rhiz_sizes))


## |-|- Create list of bins ------------------------------

ns_list = list()
log_ns = list()
for(name in names(rhiz$zoo_files)) {
  ns_list[[name]] <- calc_nbss_width(rhiz$zoo_files[[name]],
                               vol = as.numeric(vol_s$volume[vol_s$profileid == name]))
  
  log_ns[[name]] <- calc_nbss_log(rhiz$zoo_files[[name]],
                               vol = as.numeric(vol_s$volume[vol_s$profileid == name]))
}

nsdf <- ns_list |> 
  list_to_tib('profileid')

logdf <- log_ns |>
  list_to_tib('profileid')



ns_plot <- ggplot() +
  geom_point(data = nsdf,
             aes(x = mp_size_bin,
                 y = ns)
  ) +
  labs(x = 'ESD [mm]', y = '#/L    ') +
  theme_bw()

log_plot <- ggplot() +
  geom_point(data = logdf,
             aes(x = log(mp_size_bin),
                 y = ns)
  ) +
  geom_smooth(data = logdf,
              aes(x = log(mp_size_bin),
                  y = ns),
              method = 'lm',
              se = F,
              color = 'red') +
  labs(x = 'Log(ESD [mm])', y = 'Log(#/L)') +
  theme_bw()

lm(logdf$ns ~ log(logdf$mp_size_bin))

ggarrange(ns_plot, log_plot, ncol = 2, labels = 'AUTO')
ggsave('./output/Supplement/s03_NBSS-calcs.pdf', width = 8, height = 4, units = 'in')
