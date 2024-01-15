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

# |- Interpolating Function -----------------------

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
  group_by(Date, db) |> 
  summarize(par_conc = mean(par_conc)) |> 
  ungroup() |> 
  bin_format()

par_df$Depth = par_df$mp

# |- Interpolate -----------------

par_interp <- surf_interp(par_df, 'par_conc')

par_df_early <- par_df[par_df$Date <= '2019-12-31',]
par_df_late <- par_df[par_df$Date > '2019-12-31',]

par_early_interp <- par_interp[par_interp$Date <= '2019-09-30',]
par_late_interp <- par_interp[par_interp$Date >= '2020-10-24',]

par_early_plot <- ggplot() +
  geom_tile(data = par_early_interp,
            aes(x = Date, y = Depth, fill = par_conc)) +
  geom_point(data = par_df_early,
             aes(x = as.POSIXct(Date), y = Depth),
             size = 0.25, color = 'lightgrey') +
  scale_fill_viridis_c(option = 'magma',
                       begin = 0, end = 1)+
  scale_y_reverse(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position = 'none')

par_late_plot <- ggplot() +
  geom_tile(data = par_late_interp,
            aes(x = Date, y = Depth, fill = par_conc)) +
  geom_point(data = par_df_late,
             aes(x = as.POSIXct(Date), y = Depth),
             size = 0.25, color = 'lightgrey') +
  scale_fill_viridis_c(option = 'magma',
                       begin = 0, end = 1)+
  scale_y_reverse(expand = c(0,0)) +
  theme_bw()

ggsave('./output/01_environmental/par-conc-early.pdf',
       par_early_plot,
       height = 90, width = 40, units = 'mm',dpi = 600)

ggsave('./output/01_environmental/par-conc-late.pdf',
       par_late_plot,
       height = 90, width = 150, units = 'mm',dpi = 600)

###
# CTD DATA ###########
###


ctd_data <- readRDS("./data/00_ctd-data.rds") |> 
  list_to_tib('cruise_id')


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

#clean up ctd_data

ctd_data_interp <- list()
for(value in c('temp','sal','o2','RFU')) {
  ctd_data_interp[[value]] <- surf_interp(ctd_data, value)
}


# |- Plot CTD -----------------------

temp_late_interp <- ctd_data_interp$temp[ctd_data_interp$temp$Date >= '2020-10-24',]
temp_late_df <- ctd_data[ctd_data$Date >= '2020-10-24', ]
temp_plot_late <- ggplot()+
  geom_tile(data = temp_late_interp,
            aes(x = Date, y = Depth, fill = temp)) +
  geom_point(data = temp_late_df,
             aes(x = as.POSIXct(Date), y = Depth),
             size = 0.25, color = 'lightgrey') +
  scale_y_reverse(expand = c(0,0)) +
  scale_fill_viridis_c(option = 'H')
