###
# Gam plot Creation
###

rm(list = ls())


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
    theme(legend.position = 'none', axis.text = element_text(size = 8, face = 'bold'))
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

zones = c(
  `epi` = 'Epipelagic',
  `upmeso` = 'Upper Mesopelagic',
  `lomeso` = 'Lower Mesopelagic'
)

# need to make the plot into a function to avoid lazy load issue
supp_plot <- function(taxon, zone, var) {
  plot_data <- data[[taxon]][[zone]]
  
   outplot <- ggplot() +
    geom_point(aes(x = plot_data[[var]]$obs[[var]], 
                   y = plot_data[[var]]$obs$intg),
               size = 3) +
    geom_line(aes(x = plot_data[[var]]$pred[[var]],
                  y = plot_data[[var]]$pred$fit)) +
    geom_ribbon(aes(x = plot_data[[var]]$pred[[var]],
                    ymin = plot_data[[var]]$pred$low,
                    ymax = plot_data[[var]]$pred$high),
                color = 'lightgrey',
                alpha = 0.2) +
    labs(x = var, y = 'Integrated Abundance',
         subtitle = paste0(taxon, ' ', zones[zone])) +
    theme_bw() +
    theme(legend.position = 'none', 
          axis.text = element_text(size = 8, face = 'bold'),
          plot.subtitle = element_text(size = 8))
   
   return(outplot)
}


#### ERRORS HERE ####################

for(zone in names(zones)) {
  plot_list = list()
  
  for(taxon in taxa) {
    if(zone %in% names(data[[taxon]])){
      print(paste0('Making plot for ', taxon, ' in ', zone))
      for(var in names(data[[taxon]][[zone]])) {
        plot_list[[paste0(taxon, '-', var)]] <- supp_plot(taxon, zone, var)
      }
    }
  }
  
  quartz()
  ggarrange(plotlist = plot_list)
}


for(name in plot_list) {
  quartz()
  print(plot_list[[name]])
}

# # DO NOT UNCOMMENT UNLESS YOU ARE SURE YOU WANT TO DELETE ALL OTHER PLOTS!
# 
# lapply(names(data), save_plots)

