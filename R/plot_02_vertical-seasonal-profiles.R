
##
# Plots ##########
##

# |- Plot Rhizaria Seasonality ------------------------------

# |- For Total Rhizaria ---------------------------------------

ggplot(rhiz_upmeso_sum) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = "#FDDBC7", mid = "#B2182B", high = "#FDDBC7", name = "Month",
                       midpoint = 6) +
  facet_wrap(.~taxa)+
  theme_bw()


# |- By Rhizaria Group ------------------------------------

ggplot(rhiz_epi_sum[which(rhiz_epi_sum$taxa == 'Rhizaria'),]) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = "#FDDBC7", mid = "#B2182B", high = "#FDDBC7", name = "Month",
                       midpoint = 6) +
  theme_bw()

# |- Vertical Distributions -------------------------------------


vert_plotter <- function(df) {
  plot <- ggplot(data = df) +
    geom_rect(aes(xmin = min_d, xmax = max_d,
                  ymax = mean, ymin = 0))+
    geom_errorbar(aes(x = mp, ymin = mean, ymax = mean + sd)) +
    coord_flip() +
    scale_x_reverse()+ 
    theme_bw()
  return(plot)
}

total_rhiz_plot <- vert_plotter(total_avg_rhiz)

taxa_plots <- list()
for(taxa in unique(rhiz_avg_group$group)) {
  taxa_plots[[taxa]] <- vert_plotter(rhiz_avg_group[which(rhiz_avg_group$group == taxa),])
}


for(taxa in names(taxa_plots)) {
  windows()
  print(taxa_plots[[taxa]] +
          labs(subtitle = taxa))
}
