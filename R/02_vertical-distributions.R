##
# Veritical Distributions
##

rm(list = ls())
library(ggplot2)
library(dplyr)
library(lubridate)
library(EcotaxaTools)

###
# Data prep #################
###

# |- Load in data -------
rhiz <- readRDS('./data/01_rhiz-densities.rds')
meta <- readRDS('./data/00_zoop-uvp.rds')$meta

# calc meta times
meta$month <- month(meta$sampledate)
meta$year <- year(meta$sampledate)

# drop any deep casts

deep_drop <- function(cast_df) {
  if(max(cast_df$max_d) > 1200) {
    outdf <- cast_df |> 
      as_tibble() |> 
      filter(max_d < 1200)
  } else {
    outdf <- cast_df
  }
  
  return(outdf)
}

rhiz <- rhiz |> 
  lapply(deep_drop)

# |- make a all-rhizaria dataframe --------------

sum_by_cast <- function(cast_df) {
  out_df <- cast_df |> 
    group_by(db, min_d, max_d, mp) |> 
    summarise(total_rhiz = sum(conc_m3))
}

total_rhiz <- rhiz |> 
  lapply(sum_by_cast)


###
# Averages for vertical profile ################
###

total_avg_rhiz <- total_rhiz |> 
  list_to_tib('profileid') |> 
  group_by(db, min_d, max_d, mp) |> 
  summarise(mean = mean(total_rhiz),
            sd = sd(total_rhiz))


rhiz_avg_group <- rhiz |> 
  list_to_tib('profileid') |> 
  group_by(db, min_d, max_d, mp, group) |> 
  summarize(mean = mean(conc_m3),
            sd = sd(conc_m3))



####
# Integrate for Seasonality ##########
####

# |- Set up--------------------------

# write a function to make integrations based on depths and reject non-matches

integrate_range <- function(df, low, high) {
  poss_ds <- seq(low, high, 25)[-1]
  if(!(all(poss_ds %in% df$max_d))) {
    return(NA)
  } else {
    filter_df <- df |> 
      as_tibble() |> 
      filter(max_d <= high & max_d > low) |> 
      integrate_all() |> 
      intg_to_tib()
    
    return(filter_df)
  }
}

# needed to modify for the total rhiz since integrate_all doesn't work!
integrate_range_total <- function(df, low, high, ...) {
  poss_ds <- seq(low, high, 25)[-1]
  if(!(all(poss_ds %in% df$max_d))) {
    return(NA)
  } else {
    filter_df <- df |> 
      as_tibble() |> 
      filter(max_d <= high & max_d > low)
    
    
    intg_val <- trap_integrate(filter_df$max_d, filter_df$total_rhiz,
                               low, high, ...)
    
    return(as.numeric(intg_val$value))
  }
}

# For trimming later
trim_NAs <- function(rhiz_list) {
  drop_list <- rhiz_list |> 
    sapply(is.logical)
  
  return(rhiz_list[-which(drop_list)])
}


# need to find a way to order by cruise for plotting (maybe hard code)
# while they are numerically ordered, the bloom cruise throws it off
cruise_id_convert <- function(rhiz_df) {
  
  rhiz_df$cruise_id <- rhiz_df$cruise_id |> 
    as.numeric()
  
  rhiz_df$cruise_id[which(rhiz_df$cruise_id == 20379)] <- 10379.5
  rhiz_df$cruise_id <- rhiz_df$cruise_id |> 
    factor(levels = sort(unique(rhiz_df$cruise_id)))
  
  return(rhiz_df)
}

# |- Total Calcs -----------------------------------
# |-|- Epipelagic --------------------

all_epi <- total_rhiz |> 
  lapply(integrate_range_total, 0, 200) |> 
  trim_NAs() |> 
  list_to_tib('profileid')

names(all_epi) <- c('intg', 'profileid')
all_epi$intg <- as.numeric(all_epi$intg)


all_epi <- meta |> 
  select(profileid, cruise_id) |> 
  right_join(all_epi, by = "profileid")

all_epi_sum <- all_epi |> 
  group_by(cruise_id) |> 
  summarize(mean_intg = mean(intg),
            sd_intg = sd(intg)) |> 
  left_join(unique(meta[,which(names(meta) %in% c('cruise_id', 'month', 'year'))]))

all_epi_sum <- all_epi_sum |> 
  cruise_id_convert()


# |-|- Upper Meso --------------------

all_upmeso <- total_rhiz |> 
  lapply(integrate_range_total, 200, 500) |> 
  trim_NAs() |> 
  list_to_tib('profileid')

names(all_upmeso) <- c('intg', 'profileid')
all_upmeso$intg <- as.numeric(all_upmeso$intg)


all_upmeso <- meta |> 
  select(profileid, cruise_id) |> 
  right_join(all_upmeso, by = "profileid")

all_upmeso_sum <- all_upmeso |> 
  group_by(cruise_id) |> 
  summarize(mean_intg = mean(intg),
            sd_intg = sd(intg)) |> 
  left_join(unique(meta[,which(names(meta) %in% c('cruise_id', 'month', 'year'))]))

all_upmeso_sum <- all_upmeso_sum |> 
  cruise_id_convert()


# |-|- Lower Meso --------------------

all_lomeso <- total_rhiz |> 
  lapply(integrate_range_total, 500, 1000) |> 
  trim_NAs() |> 
  list_to_tib('profileid')

names(all_lomeso) <- c('intg', 'profileid')
all_lomeso$intg <- as.numeric(all_lomeso$intg)


all_lomeso <- meta |> 
  select(profileid, cruise_id) |> 
  right_join(all_lomeso, by = "profileid")

all_lomeso_sum <- all_lomeso |> 
  group_by(cruise_id) |> 
  summarize(mean_intg = mean(intg),
            sd_intg = sd(intg)) |> 
  left_join(unique(meta[,which(names(meta) %in% c('cruise_id', 'month', 'year'))]))

all_lomeso_sum <- all_lomeso_sum |> 
  cruise_id_convert()



# |- Group Calculations ---------------------------------------
# |-|- Epipelagic ------------------------------
rhiz_epi <- rhiz |> 
  lapply(integrate_range, 0, 200) |> 
  trim_NAs() |> 
  list_to_tib('profileid')

rhiz_epi <- meta |> 
  select(profileid, cruise_id) |> 
  right_join(rhiz_epi, by = "profileid")

rhiz_epi_sum <- rhiz_epi |> 
  group_by(cruise_id, taxa) |> 
  summarize(mean_intg = mean(intg),
            sd_intg = sd(intg)) |> 
  left_join(unique(meta[,which(names(meta) %in% c('cruise_id', 'month', 'year'))]))



rhiz_epi_sum <- rhiz_epi_sum |> 
  cruise_id_convert()
  

# need to calc for more vertical zones

# |-|- Upper Meso ------------------------------
rhiz_upmeso <- rhiz |> 
  lapply(integrate_range, 200, 500) |> 
  trim_NAs() |> 
  list_to_tib('profileid')

rhiz_upmeso <- meta |> 
  select(profileid, cruise_id) |> 
  right_join(rhiz_upmeso, by = "profileid")

rhiz_upmeso_sum <- rhiz_upmeso |> 
  group_by(cruise_id, taxa) |> 
  summarize(mean_intg = mean(intg),
            sd_intg = sd(intg)) |> 
  left_join(unique(meta[,which(names(meta) %in% c('cruise_id', 'month', 'year'))]))

rhiz_upmeso_sum <- rhiz_upmeso_sum |> 
  cruise_id_convert()

# |-|- Lower Meso ------------------------------
rhiz_lomeso <- rhiz |> 
  lapply(integrate_range, 500, 1000) |> 
  trim_NAs() |> 
  list_to_tib('profileid')

rhiz_lomeso <- meta |> 
  select(profileid, cruise_id) |> 
  right_join(rhiz_lomeso, by = "profileid")

rhiz_lomeso_sum <- rhiz_lomeso |> 
  group_by(cruise_id, taxa) |> 
  summarize(mean_intg = mean(intg),
            sd_intg = sd(intg)) |> 
  left_join(unique(meta[,which(names(meta) %in% c('cruise_id', 'month', 'year'))]))

rhiz_lomeso_sum <- rhiz_lomeso_sum |> 
  cruise_id_convert()

###
# Save integrated data #####
###

# Totals
saveRDS(
  list(
    epi = all_epi,
    upmeso = all_upmeso,
    lomeso = all_lomeso
  ),
  "./data/02_integrated-all-rhiz.rds"
)

# Groups
saveRDS(list(
  epi = rhiz_epi,
  upmeso = rhiz_upmeso,
  lomeso = rhiz_lomeso
), "./data/02_integrated-taxa.rds"
)


  
##
# Plots ##########
##

# |- Plot Rhizaria Seasonality ------------------------------

# |- For Total Rhizaria ---------------------------------------

ggplot(all_lomeso_sum) +
  geom_bar(aes(x = cruise_id,  y = mean_intg, fill = month),
           stat = 'identity') +
  geom_errorbar(aes(x = cruise_id, ymin = mean_intg, ymax = mean_intg + sd_intg),
                stat = 'identity') +
  scale_fill_gradient2(low = "#FDDBC7", mid = "#B2182B", high = "#FDDBC7", name = "Month",
                       midpoint = 6) +
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