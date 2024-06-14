######
# Plotting Model Fit ########
######

rm(list = ls())
graphics.off()
library(ggplot2)
library(ggpubr)

# Load the data
mods <- readRDS('./data/04_gam-intg.rds')


plot_list = list()


title_names <- c(
  'tot_epi' = 'All Rhizaria Epipelagic',
  'tot_upmeso' = 'All Rhizaria Upper Mesopelagic',
  'tot_lomeso' = 'All Rhizaria Lower Mesopelagic',
  'acantharea_epi' = 'Acantharea Epipelagic',
  'acantharea_upmeso' = 'Acantharea Upper Mesopelagic',
  'acantharea_lomeso' = 'Acantharea Lower Mesopelagic',
  'aulacanthidae_epi' = 'Aulacanthidae Epipelagic',
  'aulacanthidae_upmeso' = 'Aulacanthidae Upper Mesopelagic',
  'aulacanthidae_lomeso' = 'Aulacanthidae Lower Mesopelagic',
  'aulosphaeridae_upmeso' = 'Aulosphaeridae Upper Mesopelagic',
  'aulosphaeridae_lomeso' = 'Aulosphaeridae Lower Mesopelagic',
  'castanellidae_epi' = 'Castanellidae Epipelagic',
  'coelodendridae_upmeso' = 'Coelodendridae Upper Mesopelagic',
  'coelodendridae_lomeso' = 'Coelodendridae Lower Mesopelagic',
  'collodaria_epi' = 'Collodaria Epipelagic',
  'foraminifera_epi' = 'Foraminifera Epipelagic'
)


  
plot_list <- lapply(
  names(mods), function(name) {
    ggplot() +
      geom_point(aes(x = predict(mods[[name]]),
                     y = mods[[name]]$y),
                 size = 1, color = 'black') +
      geom_abline(intercept = 0, slope = 1, color = 'grey', lty = 2) +
      labs(x = 'Predicted', y = 'Observed',
           subtitle = title_names[name]) +
      scale_x_continuous(limits = c(
        max(c(-50, min(predict(mods[[name]]))) - 0.1*min(predict(mods[[name]]))), 
        max(predict(mods[[name]])) + 0.1*max(predict(mods[[name]]))
      ))+
      scale_y_continuous(limits = c(
        0, 
        max(mods[[name]]$y) + 0.1*max(mods[[name]]$y)
      ))+
      theme_bw() +
      theme(plot.subtitle = element_text(size = 6))
  }
)

ggsave('./output/Supplement/s04_model-fits.pdf',
    plot = ggarrange(plotlist = plot_list),
    width = 8, height = 8, units = 'in')

