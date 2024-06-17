###
# Gam plot Creation
###

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
graphics.off()

###
# Just to show particle relationships




###
# Below is for first draft rendition #########
###


## |- Plot Function -----------------------

make_gam_plot <- function(var, data) {
  plot <- ggplot() +
    geom_point(aes(x = data[[var]]$obs[[var]], y = data[[var]]$obs$intg,
                   color = data[[var]]$obs$month),
               size = 5) +
    geom_line(aes(x = data[[var]]$pred[[var]],
                  y = data[[var]]$pred$fit)) +
    geom_ribbon(aes(x = data[[var]]$pred[[var]],
                    ymin = data[[var]]$pred$low,
                    ymax = data[[var]]$pred$high),
                color = 'lightgrey',
                alpha = 0.2) +
    labs(x = var, y = 'Integrated Abundance') +
    scale_color_gradient2(low = "#96e0fa", mid = "#B2182B", high = "#96e0fa", name = "Month",
                         midpoint = 6) +
    theme_bw() +
    theme(legend.position = 'none', 
          axis.text = element_text(size = 18, face = 'bold'),
          axis.title = element_blank(),
          panel.border = element_rect(colour = 'black', fill = NA,
                                      linewidth = 2))
  return(plot)
}


##
# Make and save plots #############
##

files <- grep('04_', dir('./data/'), value = T)

data <- list()
for(file in files) {
  name <- sub("^04_(.+?)-gam\\.rds$", '\\1', file)
  data[[name]] <- readRDS(paste0('./data/',file))
}
data = data[-which(names(data) == '04_gam-intg.rds')]

save_plots <- function(name) {
  group_data <- data[[name]]
  for(zone in names(group_data)) {
    for(var in names(group_data[[zone]])) {
      cur_plot <- make_gam_plot(var, group_data[[zone]])
      
      file_loc <- paste0('./output/05_gam-partials/05_',name,'-',zone,'-gam-',var,'.pdf')
      
      ggsave(file_loc, cur_plot,
             width = 6, height = 6, 
             units = 'in')
    }
  }
}


###
# For Supplemental Plots ##########
###

taxa = names(data)
taxa = taxa[-which(taxa == 'total')]


## |- Name mapping for plots ------------------
zones = c(
  `epi` = 'Epipelagic',
  `upmeso` = 'Upper Mesopelagic',
  `lomeso` = 'Lower Mesopelagic'
)

taxa_names = c(
  `acanth` = 'Acantharea',
  `aulacanthidae` = 'Aulacanthidae',
  `aulosphaeridae` = 'Aulosphaeridae',
  `castanellidae` = 'Castanellidae',
  `coelodendridae` = 'Coelodendridae',
  `collodaria` = 'Collodaria',
  `foram` = 'Foraminifera'
)

var_names = c(
  `temp` = 'Temperature [°C]',
  `sal` = 'Salinity [ppt]',
  `pp` = 'Primary Production [mg C m^-2 d^-1]',
  `avg_mass_flux_200` = 'Average Mass Flux [mg C m^-2 d^-1]',
  `Bact_enumb` = 'Bacterial Abundance [10^8 kg^-1]',
  `RFU` = 'Relative Fluorescence Units',
  `avg_fbn_200` = 'Average Nitrogen Flux [mg N m^-2 d^-1]',
  `Si` = 'Silicate [µmol kg^-1]',
  `o2` = 'Oxygen [µmol kg^-1]',
  `par_conc` = 'Particle Concentration [L^-1]',
  `avg_fbc_200` = 'Average Carbon Flux [mg C m^-2 d^-1]'
)


# need to make the plot into a function to avoid lazy load issue
supp_plot <- function(taxon, zone, var) {
  plot_data <- data[[taxon]][[zone]]
  
   outplot <- ggplot() +
    geom_point(aes(x = plot_data[[var]]$obs[[var]], 
                   y = plot_data[[var]]$obs$intg),
               size = 1) +
    geom_line(aes(x = plot_data[[var]]$pred[[var]],
                  y = plot_data[[var]]$pred$fit)) +
    geom_ribbon(aes(x = plot_data[[var]]$pred[[var]],
                    ymin = plot_data[[var]]$pred$low,
                    ymax = plot_data[[var]]$pred$high),
                color = 'lightgrey',
                alpha = 0.2) +
    labs(x = var_names[var], y = 'Integrated Abundance [indv. m^-2]',
         subtitle = paste0(taxa_names[taxon], ' ', zones[zone])) +
    theme_bw() +
    theme(legend.position = 'none', 
          axis.text = element_text(size = 6, face = 'bold'),
          plot.subtitle = element_text(size = 6),
          axis.title = element_text(size = 6))
   
   return(outplot)
}



plot_list = list()
for(zone in names(zones)) {
  plot_list[[zone]] = list()
  
  for(taxon in taxa) {
    if(zone %in% names(data[[taxon]])){
      print(paste0('Making plot for ', taxon, ' in ', zone))
      for(var in names(data[[taxon]][[zone]])) {
        plot_list[[zone]][[paste0(taxon, '-', var)]] <- supp_plot(taxon, zone, var)
      }
    }
  }
}


for(name in names(plot_list)) {
  
  ggsave(paste0('./output/Supplement/s05_gam-partials/s05_', name, '.pdf'),
         ggarrange(plotlist = plot_list[[name]],
                   labels = 'AUTO'),
         width = 8, height = 9, units = 'in')
}

# DO NOT UNCOMMENT UNLESS YOU ARE SURE YOU WANT TO DELETE ALL OTHER PLOTS!
# 
# lapply(names(data), save_plots)

