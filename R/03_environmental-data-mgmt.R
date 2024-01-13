###
# Data Management for GAMS #
###

rm(list = ls())
library(EcotaxaTools)
library(dplyr)
library(tidyr)
library(lubridate)

##
# Data Prep ############
##

# |- Original UVP data -------------
meta <- readRDS('./data/00_zoop-uvp.rds')$meta

# |- Particle data --------------------
par_conc <- readRDS('./data/01_par-conc.rds')

par_df <- par_conc |> 
  list_to_tib('profileid')

# |- Rhizaria Data -------------------------------------------------

# bin-specific data
binned_rhiz <- readRDS('./data/01_rhiz-densities.rds')

# |-|- make a all-rhizaria dataframe --------------

sum_by_cast <- function(cast_df) {
  out_df <- cast_df |> 
    group_by(db, min_d, max_d, mp) |> 
    summarise(total_rhiz = sum(conc_m3))
}

total_rhiz <- binned_rhiz |> 
  lapply(sum_by_cast)

# integrated data

intg_total_rhiz <- readRDS('./data/02_integrated-all-rhiz.rds')
intg_taxa <- readRDS('./data/02_integrated-taxa.rds')


# |- Environmental Data --------------------------------------

# readin
bottle_data <- readRDS('./data/00_bottles.rds')
ctd_data <- readRDS('./data/00_ctd-data.rds')
flux_data <- readRDS('./data/00_flux.rds')
prod_data <- readRDS('./data/00_prod.rds')

# |-|- CTD Data ----------------------------

ctd_binner <- function(ctd_df) {
  
  #filter out too deep of observations
  ctd_df <- ctd_df |> 
    filter(depth <= 1000)
  
  ctd_df$db <- cut(ctd_df$depth,
                   seq(0,1000,25))
  
  ctd_df <- ctd_df |> 
    group_by(ctd_origfilename, db, cruise_id) |> 
    summarize(temp = mean(temp),
              sal = mean(sal),
              o2 = mean(umol_kg),
              BAC = mean(BAC),
              RFU = mean(RFU),
              PAR = mean(PAR),
    )
  
  # remove any NA rows
  if(any(is.na(ctd_df$db))) {
    ctd_df <- ctd_df |> 
      filter(!is.na(db))
  }
  ctd_df$cruise_id <- as.numeric(ctd_df$cruise_id)
  
  return(ctd_df)
}

ctd_bins <- ctd_data |> 
  lapply(ctd_binner)

# merge ctd data
ctd_df <- ctd_bins |> 
  do.call(what = rbind,)

# |-|-|- Correct the CTD Data ---------------
ctd_df$o2[which(ctd_df$o2 < 0)] <- 0
ctd_df$sal[which(ctd_df$sal < 0)] <- 0


# |-|- Bottle Data ------------------------

bottle_data$cruise_id <- bottle_data$ctd_origfilename |> 
  as.character() |> 
  substr(1,5) |> 
  as.numeric()

# replace -999 with NA

bottle_data <- bottle_data |> 
  mutate_all(~ ifelse(. == -999, NA, .)) |> 
  select_if(~!all(is.na(.)))


# remove any rows where they are all NA
drop_idx <- bottle_data[,c(18:26)] |> 
  apply(1, sum, na.rm = T)

bottle_data <- bottle_data[-which(drop_idx == 0),]

# |-|-|- Split and average within ----------------------

bot_split <- bottle_data |> 
  split(f = bottle_data$ctd_origfilename)

bottle_interpolater <- function(bot_df) {
  
  interp_names <- names(bot_df)[18:26]
  
  interp_vals <- list()
  
  for(val in interp_names) {
    drange <- seq(min(bot_df$depth), max(bot_df$depth), 1)
    
    if(all(is.na(bot_df[[val]]))) {
      interp_vals[[val]] <- NA
    } else {
      interpolator <- lin_interp(
        x = bot_df$depth,
        y = bot_df[[val]],
        min_x = min(bot_df$depth),
        max_x = max(bot_df$depth)
      )
      
      interp_vals[[val]] <-  drange |> 
        interpolator()
    }
  }
  
  out_vals <- interp_vals |> 
    data.frame()
  
  out_vals$depth <- seq(min(bot_df$depth), max(bot_df$depth), 1) |> 
    floor() #floor it for easy matching
  out_vals$ctd_origfilename <- unique(bot_df$ctd_origfilename)
  out_vals$cruise_id <- unique(bot_df$cruise_id)
  
  return(out_vals)
}

bot_vals <- bot_split |> 
  lapply(bottle_interpolater)

interpolated_avgs <- bot_vals |> 
  do.call(what = rbind,) |> 
  group_by(depth,cruise_id) |> 
  summarize(NO3 = mean(NO3, na.rm = T),
            NO2 = mean(NO2, na.rm = T),
            PO41 = mean(PO41, na.rm = T),
            Si = mean(Si, na.rm = T),
            POC = mean(POC, na.rm = T),
            PON = mean(PON, na.rm = T),
            TOC = mean(TOC, na.rm = T),
            TN = mean(TN, na.rm = T),
            Bact_enumb = mean(Bact_enumb, na.rm = T))

interpolated_avgs$db <- cut(interpolated_avgs$depth, seq(0,1000,25))

binned_bottles <- interpolated_avgs |>
  filter(!is.na(db)) |> 
  group_by(cruise_id, db) |> 
  summarize(NO3 = mean(NO3, na.rm = T),
            NO2 = mean(NO2, na.rm = T),
            PO41 = mean(PO41, na.rm = T),
            Si = mean(Si, na.rm = T),
            POC = mean(POC, na.rm = T),
            PON = mean(PON, na.rm = T),
            TOC = mean(TOC, na.rm = T),
            TN = mean(TN, na.rm = T),
            Bact_enumb = mean(Bact_enumb, na.rm = T)) 


# |-|- Flux Data -------------------
flux_predict <- flux_data |> 
  select(cruise_id, depth, avg_mass_flux, avg_fbc, avg_fbn) |> 
  pivot_wider(names_from = depth, values_from = c(avg_mass_flux, avg_fbc, avg_fbn))

# |-|- Primary Productivity ----------------
prod_split <- prod_data |> 
  split(f = prod_data$cruise_id)


# |-|-|- Interpolated Productivity --------------------

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



# |-|-|- Binned Productivity ------------------------


interp_prod$db <- cut(interp_prod$depth, breaks = seq(0,150,25))

prod_bins <- interp_prod |> 
  group_by(cruise_id, db) |> 
  summarize(pp = mean(pp))

prod_bins <- prod_bins[-which(is.na(prod_bins$db)),] |> 
  bin_format()

prod_bins$cruise_id <- as.numeric(prod_bins$cruise_id)
  
# |-|-|- Integrated Euphotic Productivity ---------------------

prod_integrate <- function(pp) {
  
  drange <- c(0:140)[-which(is.na(pp))]
  pp <- pp[-which(is.na(pp))]
  out <- trap_integrate(x = drange, y = pp, min_x = min(drange), max_x = max(drange))
  return(out$value)
}


integrated_prod <- interp_prod |> 
  group_by(cruise_id) |> 
  summarize(total_prod = prod_integrate(pp))


###
# Merge Environmental And Rhizaria #############
###

# |- Total Rhizaria -----------------------




tot_rhiz_envir <- total_rhiz |>
  list_to_tib('profileid') |>
  left_join(
    par_df, by = c('profileid','db')
  ) |>
  left_join(
    meta |>
      select(profileid, ctd_origfilename),
    by = 'profileid'
  ) |>
  left_join(ctd_df,
            by = c('ctd_origfilename', 'db')) |>
  left_join(
    binned_bottles |>
      select('db', 'cruise_id','Si','Bact_enumb', 'NO3')
  ) |>
  left_join(
    flux_predict |>
      select(cruise_id, avg_mass_flux_200)
  ) |> left_join(
    prod_bins |>
      select(cruise_id, pp, db)
  )

tot_rhiz <- list()
tot_rhiz[['epi']] <- tot_rhiz_envir[tot_rhiz_envir$max_d <= 200, ]
tot_rhiz[['meso']] <- tot_rhiz_envir[tot_rhiz_envir$max_d > 200,]

# |- Individual taxa dfs -----------------------------

taxa_rhiz <- binned_rhiz |>
  list_to_tib('profileid')  |>
  left_join(
    par_df, by = c('profileid','db')
  ) |>
  left_join(
    meta |>
      select(profileid, ctd_origfilename),
    by = 'profileid'
  ) |>
  left_join(ctd_df,
            by = c('ctd_origfilename', 'db')) |>
  left_join(
    binned_bottles |>
      select('db', 'cruise_id','Si','Bact_enumb', 'NO3')
  ) |>
  left_join(
    flux_predict |>
      select(cruise_id, avg_mass_flux_200)
  ) |> left_join(
    prod_bins |>
      select(cruise_id, pp, db)
  )

taxa_rhiz <- taxa_rhiz |> split(f = taxa_rhiz$group)

keeper_names <- c('Acantharea','Aulacanthidae','Aulosphaeridae','Castanellidae',
                  'Coelodendridae','Collodaria','Foraminifera')

taxa_data <- list()
for(name in keeper_names) {
  taxa_data[[name]] <- list()
  taxa_data[[name]][["epi"]] <- taxa_rhiz[[name]] |>
    filter(max_d <= 200)
  taxa_data[[name]][['meso']] <- taxa_rhiz[[name]] |>
    filter(max_d > 200)
}

###
# Matching Integrated Data #######
###

# |- Getting enviornmental data to match ------------
ctd_df <- ctd_df |> bin_format()
ctd_df$zone <- NULL
for(r in 1:nrow(ctd_df)) {
  if(ctd_df$max_d[r] <= 200) {
    ctd_df$zone[r] <- 'epi'
  } else if( ctd_df$max_d[r] <= 500) {
    ctd_df$zone[r] <- 'upmeso'
  } else {
    ctd_df$zone[r] <- 'lomeso'
  }
}

ctd_avg_zone <- ctd_df |> 
  group_by(ctd_origfilename, zone) |> 
  summarize(temp = mean(temp),
            sal = mean(sal),
            o2 = mean(o2),
            RFU = mean(RFU))

 # |-|- Bottle Data ----------------------------
binned_bottles <- binned_bottles |> bin_format()
binned_bottles$zone <- NULL
for(r in 1:nrow(binned_bottles)) {
  if(binned_bottles$max_d[r] <= 200) {
    binned_bottles$zone[r] <- 'epi'
  } else if( binned_bottles$max_d[r] <= 500) {
    binned_bottles$zone[r] <- 'upmeso'
  } else {
    binned_bottles$zone[r] <- 'lomeso'
  }
}

intg_bot <- binned_bottles |> 
  group_by(zone, cruise_id) |> 
  summarize(NO3 = mean(NO3),
            Si = mean(Si),
            Bact_enumb = sum(Bact_enumb))
intg_bot$cruise_id <- as.character(intg_bot$cruise_id)

# |-|- Particle Data --------------------------------
par_df$zone <- NA
par_df <- par_df |> bin_format()

for(r in 1:nrow(par_df)) {
  if(par_df$max_d[r] <= 200) {
    par_df$zone[r] <- 'epi'
  } else if( par_df$max_d[r] <= 500) {
    par_df$zone[r] <- 'upmeso'
  } else {
    par_df$zone[r] <- 'lomeso'
  }
}

part_intg <- par_df |> 
  group_by(zone, profileid) |> 
  summarize(par_conc = sum(par_conc, na.rm = T))

# |-|- Productivity Data ---------------------------

prod_intg <- prod_bins |> 
  group_by(cruise_id) |> 
  summarize(pp = sum(pp, na.rm = T))

prod_intg$cruise_id <- prod_intg$cruise_id |> as.character()

flux_predict$cruise_id <- as.character(flux_predict$cruise_id)

# |- Bring it together ------------------------

# |-|- All Rhizaria ----------------------------

meta$month <- month(meta$sampledate)

all_intg <- list() 
for(zone in names(intg_total_rhiz)) {
  all_intg[[zone]] <- intg_total_rhiz[[zone]] |>
    left_join(
      part_intg[which(part_intg$zone == zone),] |> 
        ungroup() |> 
        select(profileid, par_conc), 
      by = c('profileid')
    ) |>
    left_join(
      meta |>
        select(profileid, ctd_origfilename,  month),
      by = 'profileid'
    ) |>
    left_join(
      ctd_avg_zone[which(ctd_avg_zone$zone == zone),] |> 
        select(ctd_origfilename, temp, sal, o2, RFU),
      by = c('ctd_origfilename')
    ) |>
    left_join(
      intg_bot[which(intg_bot$zone == zone),] |>
        ungroup() |> 
        select('cruise_id','Si','Bact_enumb', 'NO3')
    ) |>
    left_join(
      flux_predict |>
        select(cruise_id, avg_mass_flux_200, avg_fbc_200, avg_fbn_200)
    ) |> 
    left_join(
      prod_intg
    )
}



intg_taxa_group <- intg_taxa |> 
  lapply(function(x) split(x, f = x$taxa))
taxa_intg <- list()
for(zone in names(intg_taxa_group)) {
  for(taxa in names(intg_taxa_group[[zone]])) {
    taxa_intg[[zone]][[taxa]] <- intg_taxa_group[[zone]][[taxa]] |>
      left_join(
        part_intg[which(part_intg$zone == zone),] |> 
          ungroup() |> 
          select(profileid, par_conc), 
        by = c('profileid')
      ) |>
      left_join(
        meta |>
          select(profileid, ctd_origfilename, month),
        by = 'profileid'
      ) |>
      left_join(
        ctd_avg_zone[which(ctd_avg_zone$zone == zone),] |> 
          select(ctd_origfilename, temp, sal, o2, RFU),
        by = c('ctd_origfilename')
      ) |>
      left_join(
        intg_bot[which(intg_bot$zone == zone),] |>
          ungroup() |> 
          select('cruise_id','Si','Bact_enumb', 'NO3')
      ) |>
      left_join(
        flux_predict |>
          select(cruise_id, avg_mass_flux_200, avg_fbc_200, avg_fbn_200)
      ) |> 
      left_join(
        prod_intg
      )
  }
}



###
# Save Data ###############
###

# |- Rhizaria Data for models ---------------------------------

# Binned data
saveRDS(tot_rhiz, './data/03_total-rhizaria.rds')
saveRDS(taxa_data, './data/03_taxa-rhizaria.rds')

# Integrated data
saveRDS(all_intg, './data/03_all-integrated.rds')
saveRDS(taxa_intg, './data/03_taxa-integrated.rds')

# |- Data for plotting environmental data -----------------

# CTD data is already saved as is.
# Particle data can be formatted later
# flux data is good to go as 00

# interpolated bottles in 1m for plotting
saveRDS(interpolated_avgs, './data/s03_interpolated-avg-bottle.rds')

# interpolated productivity data
saveRDS(interp_prod, './data/s03_interpolated-productivity.rds')
