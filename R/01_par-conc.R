###
# Calculate Particle Concentraiton ##
###

rm(list = ls())

library(EcotaxaTools)
library(dplyr)

uvp_data <- readRDS('./data/00_zoop-uvp.rds') 

par_conc <- uvp_data |> 
  uvp_par_conc(min_esd = 0.184, max_esd = 0.9)

db_cutter <- function(pdf) {
  pdf$db <- cut(pdf$depth, seq(0,1000,25))
  
  pdf <- pdf |> 
    group_by(db) |> 
    summarize(par_conc = mean(par_conc))
  
  return(pdf)
}

bin_par_conc <- par_conc |> 
  lapply(db_cutter)

saveRDS(bin_par_conc,'./data/01_par-conc.rds')
