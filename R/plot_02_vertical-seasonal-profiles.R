##
# Vertical and Seasonal plots #
##

###
# Data Prep #####
###

rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggpubr)
source('./R/utils.R')

avg_bins <- readRDS('./data/02_avg-bins.rds')
avg_intg <- readRDS("./data/02_avg-integrate.rds")


###
# Vertical Plots ######
###


vert_plotter <- function(df, color) {
  plot <- ggplot(data = df) +
    geom_rect(aes(xmin = min_d, xmax = max_d,
                  ymax = mean, ymin = 0),
              fill = color,
              color = 'black')+
    geom_errorbar(aes(x = mp, ymin = mean, ymax = mean + sd)) +
    coord_flip() +
    scale_x_reverse()+ 
    labs(y = 'Indv. /m3', x = 'Depth') +
    theme_bw()
  return(plot)
}

total_rhiz_plot <- vert_plotter(avg_bins$all, 'grey')
ggsave('./output/02_vertical-profiles/all-rhizaria.pdf',total_rhiz_plot,
       dpi = 300,
       width = 4, height = 6)

taxa_plots <- list()
for(taxa in unique(avg_bins$taxa$group)) {
  cur_plot <- vert_plotter(avg_bins$taxa[which(avg_bins$taxa$group == taxa),],
                                     color = taxa_colors[taxa])
  
  ggsave(paste0('./output/02_vertical-profiles/',taxa,'.pdf'),
         cur_plot, dpi = 600, width = 4, height = 6)
}


####
# Seasonality ##########
####

# |- For All Rhizaria -------------------
all_epi <- ggplot(avg_intg$all$epi) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = seasonal_scale['winter'],
                       mid = seasonal_scale['summer'],
                       high = seasonal_scale['winter'],
                       name = "Month",
                       midpoint = 6) +
  theme_bw()

ggsave('./output/02_seasonality/all-epi.pdf',
       all_epi,
       width = 190, height = 71.25, dpi = 600, units = 'mm')

all_upmeso <- ggplot(avg_intg$all$upmeso) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = seasonal_scale['winter'],
                       mid = seasonal_scale['summer'],
                       high = seasonal_scale['winter'],
                       name = "Month",
                       midpoint = 6) +
  theme_bw()

ggsave('./output/02_seasonality/all-upmeso.pdf',
       all_upmeso,
       width = 190, height = 71.25, dpi = 600, units = 'mm')


all_lomeso <- ggplot(avg_intg$all$lomeso) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = seasonal_scale['winter'],
                       mid = seasonal_scale['summer'],
                       high = seasonal_scale['winter'],
                       name = "Month",
                       midpoint = 6) +
  theme_bw()

ggsave('./output/02_seasonality/all-lomeso.pdf',
       all_lomeso,
       width = 190, height = 71.25, dpi = 600, units = 'mm')

# |- For Taxa-Specific Rhizaria ---------------------------------------

epi_sub1 <- c('Acantharea','Collodaria','Castanellidae', 
              'Foraminifera','Aulacanthidae','Aulosphaeridae')

epi_taxa <- ggplot(avg_intg$taxa$epi[which(avg_intg$taxa$epi$taxa %in% epi_sub1),]) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = seasonal_scale['winter'],
                       mid = seasonal_scale['summer'],
                       high = seasonal_scale['winter'],
                       name = "Month",
                       midpoint = 6) +
  facet_wrap(.~taxa, scales = 'free_y')+
  theme_bw()

ggsave('./output/02_seasonality/epi_taxa.pdf',
       epi_taxa, width = 190, height = 120, units = 'mm',
       dpi = 600)

upmeso_sub <- c('Acantharea','Aulacanthidae','Aulosphaeridae',
                'Coelodendridae','Collodaria','Foraminifera',
                'Rhizaria')

upmeso_taxa <- ggplot(avg_intg$taxa$upmeso[which(avg_intg$taxa$upmeso$taxa %in% upmeso_sub),]) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = seasonal_scale['winter'],
                       mid = seasonal_scale['summer'],
                       high = seasonal_scale['winter'],
                       name = "Month",
                       midpoint = 6) +
  facet_wrap(.~taxa, scales = 'free_y')+
  theme_bw()

ggsave('./output/02_seasonality/upmeso_taxa.pdf',
       upmeso_taxa, width = 190, height = 190, units = 'mm',
       dpi = 600)


lomeso_sub <- c('Acantharea','Aulacanthidae','Aulosphaeridae',
                'Coelodendridae','Foraminifera',
                'Rhizaria')

lomeso_taxa <- ggplot(avg_intg$taxa$lomeso[which(avg_intg$taxa$lomeso$taxa %in% lomeso_sub),]) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = seasonal_scale['winter'],
                       mid = seasonal_scale['summer'],
                       high = seasonal_scale['winter'],
                       name = "Month",
                       midpoint = 6) +
  facet_wrap(.~taxa, scales = 'free_y')+
  theme_bw()

ggsave('./output/02_seasonality/lomeso_taxa.pdf',
       lomeso_taxa, width = 190, height = 120, units = 'mm',
       dpi = 600)

# |- Stacked Contribution --------------------
epi_contrib <- ggplot(avg_intg$taxa$epi) +
  geom_bar(aes(x = cruise_id, y = mean_intg,
               fill = taxa),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values = taxa_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_blank())

upmeso_contrib <- ggplot(avg_intg$taxa$upmeso) +
  geom_bar(aes(x = cruise_id, y = mean_intg,
               fill = taxa),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values = taxa_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_blank())

lomeso_contrib <- ggplot(avg_intg$taxa$lomeso) +
  geom_bar(aes(x = cruise_id, y = mean_intg,
               fill = taxa),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values = taxa_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_blank())


ggsave('./output/02_seasonality/epi-contrib.pdf',
       epi_contrib,
       width = 75, height = 26, dpi = 500, units = 'mm')
ggsave('./output/02_seasonality/upmeso-contrib.pdf',
       upmeso_contrib,
       width = 75, height = 26, dpi = 500, units = 'mm')
ggsave('./output/02_seasonality/lomeso-contrib.pdf',
       lomeso_contrib,
       width = 75, height = 26, dpi = 500, units = 'mm')
