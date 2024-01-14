###
# Gam plot Creation
###

rm(list = ls())




# |- Plot Function -----------------------

make_gam_plot <- function(var, data) {
  plot <- ggplot() +
    geom_point(aes(x = data[[var]]$obs[[var]], y = data[[var]]$obs$intg,
                   color = data[[var]]$obs$month)) +
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

save_plots <- function(name) {
  group_data <- data[[name]]
  for(zone in names(group_data)) {
    for(var in names(group_data[[zone]])) {
      cur_plot <- make_gam_plot(var, group_data[[zone]])
      
      file_loc <- paste0('./output/05_',name,'-',zone,'-gam-',var,'.pdf')
      
      ggsave(file_loc, cur_plot,
             width = 6, height = 6, 
             units = 'in')
    }
  }
}

# # DO NOT UNCOMMENT UNLESS YOU ARE SURE YOU WANT TO DELETE ALL OTHER PLOTS!
# 
# lapply(names(data), save_plots)

