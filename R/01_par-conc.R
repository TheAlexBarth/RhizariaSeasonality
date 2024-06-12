###
# Calculate Particle Concentraiton ##
###

rm(list = ls())

library(EcotaxaTools)
library(dplyr)

uvp_data <- readRDS('./data/00_zoop-uvp.rds') 

par_conc <- uvp_data |> 
  uvp_par_conc(bin_limits = c(0.184, 0.450, 0.900))

# rename particles sizes
renamer <- function(x) {
  switch(x,
         '(0.184,0.45]' = 'small',
         '(0.45,0.9]' = 'large')
}

# loop through and apply
for(cast in names(par_conc)) {
  par_conc[[cast]][['esd_bin']] <- par_conc[[cast]][['esd_bin']]|> 
    sapply(renamer)
}


db_cutter <- function(pdf) {
  pdf$db <- cut(pdf$depth, seq(0,1000,25))
  
  pdf <- pdf |> 
    group_by(db, esd_bin) |> 
    summarize(par_conc = mean(par_conc))
  
  return(pdf)
}

bin_par_conc <- par_conc |> 
  lapply(db_cutter)

saveRDS(bin_par_conc,'./data/01_par-conc.rds')
