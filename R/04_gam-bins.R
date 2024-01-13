####
# Running the GAMS ###
####
rm(list = ls())

library(mgcv)
library(EcotaxaTools)
library(dplyr)
library(tidyr)

##
# Load in Data #########
##

tot_rhiz <- readRDS('./data/03_total-rhizaria.rds')
taxa_rhiz <- readRDS('./data/03_taxa-rhizaria.rds')

###
# Total Rhizaria Model ##########
###


# |- Epipelagic ----------------------
tot_mod_epi_full <- gam(
  total_rhiz ~ s(temp, k = 5) + s(sal, k = 5) + s(o2, k = 5) + s(par_conc, k = 5) +
    s(RFU, k = 5) + s(Bact_enumb, k = 5) + s(pp, k = 5) + s(avg_mass_flux_200, k = 5) +
    s(Si, k = 5),
  data = tot_rhiz$epi, method = 'ML', select = T
)
tot_mod_epi_full |> summary()

#backwards selection was done iteratively over several reduced models
tot_mod_epi_red <- gam(
  total_rhiz ~ s(temp, k = 5) + s(sal, k = 5) + s(par_conc, k = 5) + 
    s(Bact_enumb, k = 5) + s(avg_mass_flux_200, k = 5),
  data = tot_rhiz$epi, method = 'ML', select = T
)
tot_mod_epi_red |> summary()

# |- Mesopelagic -----------------------------------

tot_mod_meso_full <- gam(
  total_rhiz ~  s(temp, k = 5) + s(sal, k = 5) + s(o2, k = 5) + s(par_conc, k = 5) +
    s(Bact_enumb, k = 5) + s(avg_mass_flux_200, k = 5) +
    s(Si, k = 5),
  data = tot_rhiz$meso, method = 'ML', select = T
)
tot_mod_meso_full |> summary()

#backwards selection was done iteratively over several reduced models
tot_mod_meso_red <- gam(
  total_rhiz ~   s(sal, k = 5) + s(o2, k = 5) + s(par_conc, k = 5) +
    s(Bact_enumb, k = 5) + s(avg_mass_flux_200, k = 5) +
    s(Si, k = 5),
  data = tot_rhiz$meso, method = 'ML', select = T
)
tot_mod_meso_red |> summary()



###
# Taxa Specific Models ################
###

# |- Acantharea ---------------------------

acanth_epi_full <- gam(
  conc_m3 ~ s(temp, k = 5) + s(sal, k = 5) + s(o2, k = 5) + s(RFU, k = 5) +
    s(par_conc, k = 5) + s(pp, k = 5) + s(Bact_enumb, k = 5) + s(avg_mass_flux_200, k=5),
  data = taxa_rhiz$Acantharea$epi, method = 'ML', select = T
)
acanth_epi_full |> summary()


acanth_epi_red <- gam(
  conc_m3 ~ s(temp, k = 5) + s(sal, k = 5) + s(o2, k = 5) +
    s(par_conc, k = 5) + s(Bact_enumb, k = 5) + s(avg_mass_flux_200, k=5),
  data = taxa_rhiz$Acantharea$epi, method = 'ML', select = T
)
acanth_epi_red |> summary()

acanth_meso_full <- gam(
  conc_m3 ~ s(temp, k = 5) + s(sal, k = 5) + s(o2, k = 5) +
    s(par_conc, k = 5) + s(Bact_enumb, k = 5) + s(avg_mass_flux_200, k=5),
  data = taxa_rhiz$Acantharea$meso, method = 'ML', select = T
)
acanth_meso_full |> summary()

acanth_meso_red <- gam(
  conc_m3 ~ s(temp, k = 5) + s(sal, k = 5) + s(o2, k = 5) +
    s(par_conc, k = 5) + s(avg_mass_flux_200, k=5),
  data = taxa_rhiz$Acantharea$meso, method = 'ML', select = T
)
acanth_meso_red |> summary()
  

# |- Aulacanthidae --------------------------------------
aulacan_epi <- gam(
  conc_m3 ~ s(temp, k = 5) + s(sal, k = 5) +
    s(avg_mass_flux_200, k=5),
  data = taxa_rhiz$Aulacanthidae$epi, method = 'ML', select = T
)
aulacan_epi |> summary()


aulacan_meso <- gam(
  conc_m3 ~ s(temp, k = 5)+ s(par_conc, k = 5) +
    s(avg_mass_flux_200, k = 5) +
    s(Si, k = 5),
  data = taxa_rhiz$Aulacanthidae$meso, method = 'ML', select = T
)
aulacan_meso |>  summary()

# |- Aulosphaeridae -----------------------------------------------
