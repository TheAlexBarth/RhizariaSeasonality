###
# Plotting Environemntal Output
###

rm(list = ls())
library(EcotaxaTools)
library(MBA)
library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(scales)
source('./R/utils.R')

## |- Interpolating Function -----------------------

surf_interp <- function(df,
                        value,
                        res = 300,
                        ...) {
  xyz <- cbind(decimal_date(df$Date),
               df$Depth,
               df[[value]])
  interp_output <- mba.surf(xyz, no.X = res, no.Y = res, ...)
  long_interp <- interp_output$xyz.est$z |> 
    as.data.frame() |> 
    pivot_longer(cols = everything())
  
  outdf <- data.frame(Date = date_decimal(interp_output$xyz.est$x),
                      Depth = interp_output$xyz.est$y) |> 
    expand.grid()
  outdf[[value]] <- c(interp_output$xyz.est$z)
  return(outdf)
}

##
# Particle Data #####
##
uvp_meta <- readRDS('./data/00_zoop-uvp.rds')$meta
par_conc <- readRDS('./data/01_par-conc.rds')

par_df <- par_conc |> 
  list_to_tib('profileid') |> 
  left_join(
    uvp_meta |> 
      select(profileid,sampledate)
  )

par_df$Date <- as.Date(par_df$sampledate)

# for casts occuring on the same date, mean
par_df <- par_df |> 
  group_by(Date, db, esd_bin) |> 
  summarize(par_conc = mean(par_conc)) |> 
  ungroup() |> 
  bin_format()

par_df$Depth = par_df$mp

## |- Interpolate -----------------

## |-|- NEEDS REVISION FOR SMALL & LARGE --------------------------

## |-|- Small -----------------------------------
par_df <- par_df |> 
  pivot_wider(names_from = 'esd_bin',
              values_from = 'par_conc')


small_par_interp <- surf_interp(par_df, 'small')
small_par_interp$small[small_par_interp$small > 1] = 1

par_df_early <- par_df[par_df$Date <= '2019-12-31',]
par_df_late <- par_df[par_df$Date > '2019-12-31',]

small_par_early_interp <- small_par_interp[small_par_interp$Date <= '2019-09-30',]
small_par_late_interp <- small_par_interp[small_par_interp$Date >= '2020-10-24',]

small_par_early_plot <- ggplot() +
  geom_tile(data = small_par_early_interp,
            aes(x = Date, y = Depth, fill = small)) +
  geom_point(data = par_df_early,
             aes(x = as.POSIXct(Date), y = Depth),
             size = 0.1, color = 'lightgrey', alpha = 0.75) +
  scale_fill_viridis_c(option = 'magma',
                       begin = 0, end = 1)+
  scale_y_reverse(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = 2))

small_par_late_plot <- ggplot() +
  geom_tile(data = small_par_late_interp,
            aes(x = Date, y = Depth, fill = small)) +
  geom_point(data = par_df_late,
             aes(x = as.POSIXct(Date), y = Depth),
             size = 0.1, color = 'lightgrey', alpha = 0.75) +
  scale_fill_viridis_c(option = 'magma',
                       begin = 0, end = 1,
                       breaks = c(0.2, 0.4, 0.6))+
  scale_y_reverse(expand = c(0,0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 2))

ggsave('./output/01_environmental/small-par-conc-early.pdf',
       small_par_early_plot,
       height = 40, width = 20, units = 'mm',dpi = 600)

ggsave('./output/01_environmental/small-par-conc-late.pdf',
       small_par_late_plot,
       height = 40, width = 70, units = 'mm',dpi = 600)


## |-|- Large ----------------------------------------


large_par_interp <- surf_interp(par_df, 'large')

large_par_interp$large[large_par_interp$large > 1] = 1
large_par_early_interp <- large_par_interp[large_par_interp$Date <= '2019-09-30',]
large_par_late_interp <- large_par_interp[large_par_interp$Date >= '2020-10-24',]

large_par_early_plot <- ggplot() +
  geom_tile(data = large_par_early_interp,
            aes(x = Date, y = Depth, fill = large)) +
  geom_point(data = par_df_early,
             aes(x = as.POSIXct(Date), y = Depth),
             size = 0.1, color = 'lightgrey', alpha = 0.75) +
  scale_fill_viridis_c(option = 'magma',
                       begin = 0, end = 1)+
  scale_y_reverse(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = 2))

large_par_late_plot <- ggplot() +
  geom_tile(data = large_par_late_interp,
            aes(x = Date, y = Depth, fill = large)) +
  geom_point(data = par_df_late,
             aes(x = as.POSIXct(Date), y = Depth),
             size = 0.1, color = 'lightgrey', alpha = 0.75) +
  scale_fill_viridis_c(option = 'magma',
                       begin = 0, end = 1,
                       breaks = c(0.2,0.4,0.6))+
  scale_y_reverse(expand = c(0,0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 2))

ggsave('./output/01_environmental/large-par-conc-early.pdf',
       large_par_early_plot,
       height = 40, width = 20, units = 'mm',dpi = 600)

ggsave('./output/01_environmental/large-par-conc-late.pdf',
       large_par_late_plot,
       height = 40, width = 70, units = 'mm',dpi = 600)


###
# CTD DATA ###########
###


ctd_data <- readRDS("./data/00_ctd-data.rds") |> 
  list_to_tib('cruise_id')

## |- Prep -----------------------------------

ctd_data <- ctd_data[which(ctd_data$ctd_origfilename %in% uvp_meta$ctd_origfilename),]

ctd_data <- ctd_data[-which(ctd_data$depth == -999),]
ctd_data <- ctd_data[-which(ctd_data$RFU == -999),]
ctd_data <- ctd_data[-which(ctd_data$umol_kg == -999),]

ctd_data <- ctd_data |> 
  filter(depth <=1000) |> 
  group_by(date(datetime), depth) |>
  summarize(temp = mean(temp, na.rm = T),
            sal = mean(sal, na.rm = T),
            o2 = mean(umol_kg, na.rm = T),
            RFU = mean(RFU, na.rm = T)) |> 
  ungroup()

ctd_data$Depth <- ctd_data$depth
ctd_data$Date <- ctd_data$`date(datetime)`

## |-|- Interpolate ctd ----------------
ctd_data_interp <- list()
for(value in c('temp','sal','o2','RFU')) {
  ctd_data_interp[[value]] <- surf_interp(ctd_data, value)
}


## |- Plot CTD -----------------------

ctd_plotter <- function(value, color_option) {
  
  
  
  late_interp <- ctd_data_interp[[value]][ctd_data_interp[[value]]$Date >= '2020-10-24',]
  late_df <- ctd_data[ctd_data$Date >= '2020-10-24', ]
  
  ctd_plot_late <- ggplot()+
    geom_tile(data = late_interp,
              aes(x = Date, y = Depth, fill = late_interp[[value]])) +
    geom_point(data = late_df[seq(1,nrow(late_df),10),],
               aes(x = as.POSIXct(Date), y = Depth),
               size = 0.005, color = 'lightgrey', alpha = 0.75) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_viridis_c(option = color_option) +
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 2))
  
  
  early_interp <- ctd_data_interp[[value]][ctd_data_interp[[value]]$Date <= '2019-09-30',]
  early_df <- ctd_data[ctd_data$Date <= '2019-09-30',]
  
  ctd_plot_early <- ggplot()+
    geom_tile(data = early_interp,
              aes(x = Date, y = Depth, fill = early_interp[[value]])) +
    geom_point(data = early_df[seq(1,nrow(early_df),10),],
               aes(x = as.POSIXct(Date), y = Depth),
               size = 0.005, color = 'lightgrey', alpha = 0.75) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_viridis_c(option = color_option) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          axis.text = element_text(size = 2))
  
  early_title <- paste0('./output/01_environmental/', value, '-early.pdf')
  late_title <- paste0('./output/01_environmental/', value, '-late.pdf')
  ggsave(early_title,
         ctd_plot_early,
         height = 40, width = 20, units = 'mm',dpi = 600)
  
  ggsave(late_title,
         ctd_plot_late,
         height = 40, width = 70, units = 'mm',dpi = 600)
  
}

for(value in c('temp','sal', 'o2', 'RFU')) {
  color_option <- switch(value,
                         "temp" = 'H',
                         'sal' = 'E',
                         'o2' = 'G',
                         'RFU' = 'D')
  ctd_plotter(value, color_option)
}


###
# Bottle Data ############
##

bot_data <- readRDS('./data/s03_interpolated-avg-bottle.rds')
raw_bot <- readRDS('./data/00_bottles.rds')

bot_data <- bot_data |> 
  filter(depth <= 1000)

raw_bot <- raw_bot |> 
  filter(depth <= 1000)

## |- Prep -------------------
uvp_cruise_meta <- uvp_meta |> 
  group_by(cruise_id = as.numeric(cruise_id)) |> 
  summarize(Date = mean(sampledate))

bot_data <- bot_data |> 
  select(Depth = depth,
         cruise_id, NO3, Si, Bact = Bact_enumb) |> 
  left_join(uvp_cruise_meta) |> 
  filter(!is.nan(Si),
         !is.nan(NO3),
         !is.nan(Bact))


raw_bot$cruise_id <- raw_bot$ctd_origfilename |> 
  as.character() |> 
  substr(1,5) |> 
  as.numeric()

  
raw_bot <- raw_bot |> 
  left_join(uvp_cruise_meta) |> 
  select(Depth = depth,
         NO3, Si,
         Bact = Bact_enumb,
         cruise_id, Date)

raw_bot_list <- list()
bot_vals <- c("NO3",'Si','Bact')
for(val in bot_vals) {
  raw_bot_list[[val]] <- raw_bot |> 
    select(Depth, val, cruise_id, Date)
  
  raw_bot_list[[val]] <- raw_bot_list[[val]][which(raw_bot_list[[val]][[val]] != -999),]
}



## |-|- Interpolate Bottle data -------------------------

bot_interp <- list()
for(value in bot_vals) {
  bot_interp[[value]] <- surf_interp(bot_data, value)
}


## |- Bot Plots ---------------------

bot_plotter <- function(value, color_option) {
  
  late_interp <- bot_interp[[value]][bot_interp[[value]]$Date >= '2020-10-24',]
  late_df <- raw_bot_list[[value]][raw_bot_list[[value]]$Date >= '2020-10-24', ]
  
  bot_plot_late <- ggplot()+
    geom_tile(data = late_interp,
              aes(x = Date, y = Depth, fill = late_interp[[value]])) +
    geom_point(data = late_df,
               aes(x = as.POSIXct(Date), y = Depth),
               size = 0.25, color = 'lightgrey', alpha = 0.75) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_viridis_c(option = color_option) +
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 2))
  
  early_interp <- bot_interp[[value]][bot_interp[[value]]$Date <= '2019-09-30',]
  early_df <- raw_bot_list[[value]][raw_bot_list[[value]]$Date <= '2019-09-30',]
  
  bot_plot_early <- ggplot()+
    geom_tile(data = early_interp,
              aes(x = Date, y = Depth, fill = early_interp[[value]])) +
    geom_point(data = early_df,
               aes(x = as.POSIXct(Date), y = Depth),
               size = 0.25, color = 'lightgrey', alpha = 0.75) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_viridis_c(option = color_option) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          axis.text = element_text(size = 2))
  
  early_title <- paste0('./output/01_environmental/', value, '-early.pdf')
  late_title <- paste0('./output/01_environmental/', value, '-late.pdf')
  ggsave(early_title,
         bot_plot_early,
         height = 40, width = 20, units = 'mm',dpi = 600)
  
  ggsave(late_title,
         bot_plot_late,
         height = 40, width = 70, units = 'mm',dpi = 600)
}

for(value in bot_vals) {
  color_option <- switch(value,
                         'NO3' = 'B',
                         'Si' = 'C',
                         "Bact" = 'F')
  
  bot_plotter(value, color_option)
}


###
# Flux & PP ######
###

## |- Flux -----------------------------------

flux_data <- readRDS('./data/00_flux.rds')
flux_data <- flux_data |> 
  select(cruise_id, Depth = depth, avg_mass_flux, avg_fbc, avg_fbn) |> 
  filter(Depth == 200) |> 
  left_join(
    uvp_cruise_meta
  )

flux_data$month <- month(flux_data$Date)

flux_plotter <- function(value) {
  plot <- ggplot(flux_data) +
    geom_bar(aes(x = as.character(cruise_id), y = .data[[value]],
                 fill = month),
             stat = 'identity') +
    scale_fill_gradient2(low = seasonal_scale['winter'], mid = seasonal_scale['summer'],
                         high = seasonal_scale['winter'],
                         name = 'Month', midpoint = 6) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw()
}

flux_plots <- list()
for(value in c('avg_mass_flux','avg_fbc','avg_fbn')) {
  flux_plots[[value]] <- flux_plotter(value)
  
  ggsave(paste0('./output/01_environmental/',value, '.pdf'),
         flux_plots[[value]],
         width = 140, dpi = 500, units = 'mm')
}


## |- Productivity ----------

prod_data <- readRDS('./data/00_prod.rds')

## |-|- Primary Productivity ----------------
prod_split <- prod_data |> 
  split(f = prod_data$cruise_id)


## |-|-|- Interpolated Productivity --------------------

prod_interpolator <- function(prod_df) {
  drange = seq(0,140,1)
  
  interpolator <- lin_interp(prod_df$depth, prod_df$pp,
                             min(prod_df$depth), max(prod_df$depth))
  
  interp_pp <- interpolator(drange)
  return(
    data.frame(
      depth = drange,
      pp = interp_pp
    )
  )
}

interp_prod <- prod_split |> 
  lapply(prod_interpolator) |> 
  list_to_tib('cruise_id')



## |-|-|- Binned Productivity ------------------------


interp_prod$db <- cut(interp_prod$depth, breaks = seq(0,150,25))

prod_bins <- interp_prod |> 
  group_by(cruise_id, db) |> 
  summarize(pp = mean(pp))

prod_bins <- prod_bins[-which(is.na(prod_bins$db)),] |> 
  bin_format()

prod_bins$cruise_id <- as.numeric(prod_bins$cruise_id)

prod_intg <- prod_bins |> 
  group_by(cruise_id) |> 
  summarize(pp = sum(pp, na.rm = T))

prod_intg$cruise_id <- prod_intg$cruise_id |> as.character()


prod_intg$cruise_id <- prod_intg$cruise_id |> 
  as.numeric()



prod_intg <- prod_intg |> 
  left_join(
    uvp_cruise_meta
  )
prod_intg$month <- month(prod_intg$Date)

prod_intg$cruise_id[which(prod_intg$cruise_id == 20379)] <- 10379.5

prod_intg$cruise_id <- prod_intg$cruise_id |> 
  factor(levels = sort(unique(prod_intg$cruise_id))) |> 
  as.character()

prod_plot <- ggplot(prod_intg) +
  geom_bar(aes(x = cruise_id, y = pp, fill = month),
           stat = 'identity') +
  scale_fill_gradient2(low = seasonal_scale['winter'], mid = seasonal_scale['summer'],
                       high = seasonal_scale['winter'],
                       name = 'Month', midpoint = 6) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

ggsave('./output/01_environmental/prod.pdf',
       prod_plot,
       width = 140, dpi = 500, units =)
