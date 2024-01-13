###
# Gam on integrated data ###
###

rm(list = ls())

library(mgcv)
library(EcotaxaTools)
library(dplyr)
library(tidyr)
library(mgcViz)

##
# Load in Data #########
##

tot_rhiz <- readRDS('../data/03_all-integrated.rds')
taxa_rhiz <- readRDS('../data/03_taxa-integrated.rds')

###
# Total Rhizaria ##########
###

# |- Epi pelagic -----------------------

# full model 
full_epi_total <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = tot_rhiz$epi, method = 'ML', select = T
)
full_epi_total |>  summary()

# reduced model
red_epi_total <- gam(
  intg ~ s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) + s(NO3, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = tot_rhiz$epi, method = 'ML', select = T
)
red_epi_total |>  summary()
AIC(red_epi_total, full_epi_total)
BIC(red_epi_total, full_epi_total)

# |- UpMeso ----------------------------------

full_upmeso_tot <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) +  s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = tot_rhiz$upmeso, method = 'ML', select = T 
)
full_upmeso_tot |> summary()


red_upmeso_tot <- gam(
  intg ~
    s(Si, k = 6) + 
    s(avg_mass_flux_200, k = 6) +
    s(par_conc, k = 6),
  data = tot_rhiz$upmeso, method = 'ML', select = T 
)
red_upmeso_tot |> summary()
AIC(full_upmeso_tot, red_upmeso_tot)
BIC(full_upmeso_tot, red_upmeso_tot)

# |- LoMeso -----------------------------


full_lomeso_tot <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = tot_rhiz$lomeso, method = 'ML', select = T 
)
full_lomeso_tot |> summary()


red_lomeso_tot <- gam(
  intg ~ 
    s(Si, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = tot_rhiz$lomeso, method = 'ML', select = T 
)
red_lomeso_tot |> summary()
AIC(full_lomeso_tot, red_lomeso_tot)
BIC(full_lomeso_tot, red_lomeso_tot)


##
# Acantharea #########
##

# |- Epipelagic -----------------------

# full model 
full_epi_acantharea <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Bact_enumb, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$epi$Acantharea, method = 'ML', select = T
)
full_epi_acantharea |>  summary()

# reduced model
red_epi_acantharea <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Bact_enumb, k = 6) +
    s(avg_fbc_200, k = 6),
  data = taxa_rhiz$epi$Acantharea, method = 'ML', select = T
)
red_epi_acantharea |>  summary()
AIC(red_epi_acantharea, full_epi_acantharea)
BIC(red_epi_acantharea, full_epi_acantharea)

# |- UpMeso ----------------------------------

full_upmeso_acantharea <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Acantharea, method = 'ML', select = T 
)
full_upmeso_acantharea |> summary()


red_upmeso_acantharea <- gam(
  intg ~
    s(o2, k = 6) +
    s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Acantharea, method = 'ML', select = T 
)
red_upmeso_acantharea |> summary()
AIC(full_upmeso_acantharea, red_upmeso_acantharea)
BIC(full_upmeso_acantharea, red_upmeso_acantharea)

# |- LoMeso -----------------------------

full_lomeso_acantharea <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Acantharea, method = 'ML', select = T 
)
full_lomeso_acantharea |> summary()

red_lomeso_acantharea <- gam(
  intg ~ s(temp, k = 6) + 
    s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Acantharea, method = 'ML', select = T 
)
red_lomeso_acantharea |> summary()
AIC(full_lomeso_acantharea, red_lomeso_acantharea)
BIC(full_lomeso_acantharea, red_lomeso_acantharea)


##
# Aulacanthidae #########
##

# |- Epipelagic -----------------------

# full model 
full_epi_Aulacanthidae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$epi$Aulacanthidae, method = 'ML', select = T
)
full_epi_Aulacanthidae |>  summary()

# reduced model
red_epi_Aulacanthidae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(RFU, k = 6) +
    s(Bact_enumb, k = 6) + 
    s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$epi$Aulacanthidae, method = 'ML', select = T
)
red_epi_Aulacanthidae |>  summary()
AIC(red_epi_Aulacanthidae, full_epi_Aulacanthidae)
BIC(red_epi_Aulacanthidae, full_epi_Aulacanthidae)

# |- UpMeso ----------------------------------

full_upmeso_Aulacanthidae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Aulacanthidae, method = 'ML', select = T 
)
full_upmeso_Aulacanthidae |> summary()


red_upmeso_Aulacanthidae <- gam(
  intg ~ 
    s(NO3, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Aulacanthidae, method = 'ML', select = T 
)
red_upmeso_Aulacanthidae |> summary()
AIC(full_upmeso_Aulacanthidae, red_upmeso_Aulacanthidae)
BIC(full_upmeso_Aulacanthidae, red_upmeso_Aulacanthidae)

# |- LoMeso -----------------------------

full_lomeso_Aulacanthidae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Aulacanthidae, method = 'ML', select = T 
)
full_lomeso_Aulacanthidae |> summary()

red_lomeso_Aulacanthidae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + 
    s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Aulacanthidae, method = 'ML', select = T 
)
red_lomeso_Aulacanthidae |> summary()
AIC(full_lomeso_Aulacanthidae, red_lomeso_Aulacanthidae)
BIC(full_lomeso_Aulacanthidae, red_lomeso_Aulacanthidae)


##
# Aulosphaeridae #########
##

# |- Epipelagic -----------------------

# # full model 
# full_epi_Aulosphaeridae <- gam(
#   intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
#     s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
#     s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
#     s(pp, k = 6) + s(par_conc, k = 6),
#   data = taxa_rhiz$epi$Aulosphaeridae, method = 'ML', select = T
# )
# full_epi_Aulosphaeridae |>  summary()
# 
# # reduced model
# red_epi_Aulosphaeridae <- gam(
#   intg ~ s(sal, k = 6) + s(o2, k = 6) + 
#     s(NO3, k = 6) +
#     s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6),
#   data = taxa_rhiz$epi$Aulosphaeridae, method = 'ML', select = T
# )
# red_epi_Aulosphaeridae |>  summary()


# |- UpMeso ----------------------------------

full_upmeso_Aulosphaeridae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Aulosphaeridae, method = 'ML', select = T 
)
full_upmeso_Aulosphaeridae |> summary()


red_upmeso_Aulosphaeridae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) +
    s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Aulosphaeridae, method = 'ML', select = T 
)
red_upmeso_Aulosphaeridae |> summary()
AIC(full_upmeso_Aulosphaeridae, red_upmeso_Aulosphaeridae)
BIC(full_upmeso_Aulosphaeridae, red_upmeso_Aulosphaeridae)

# |- LoMeso -----------------------------

full_lomeso_Aulosphaeridae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Aulosphaeridae, method = 'ML', select = T 
)
full_lomeso_Aulosphaeridae |> summary()

red_lomeso_Aulosphaeridae <- gam(
  intg ~ 
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Aulosphaeridae, method = 'ML', select = T 
)
red_lomeso_Aulosphaeridae |> summary()
AIC(full_lomeso_Aulosphaeridae, red_lomeso_Aulosphaeridae)
BIC(full_lomeso_Aulosphaeridae, red_lomeso_Aulosphaeridae)

###
# Castanellidae #####
###

# full model
full_epi_Castanellidae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$epi$Castanellidae, method = 'ML', select = T
)
full_epi_Castanellidae |>  summary()

# reduced model
red_epi_Castanellidae <- gam(
  intg ~ s(temp, k = 6) + 
    s(Si, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$epi$Castanellidae, method = 'ML', select = T
)
red_epi_Castanellidae |>  summary()

AIC(full_epi_Castanellidae, red_epi_Castanellidae)
BIC(full_epi_Castanellidae, red_epi_Castanellidae)

###
# Coelodendridae ###########
###


# |- UpMeso ----------------------------------

full_upmeso_Coelodendridae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Coelodendridae, method = 'ML', select = T 
)
full_upmeso_Coelodendridae |> summary()


red_upmeso_Coelodendridae <- gam(
  intg ~  
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Coelodendridae, method = 'ML', select = T 
)
red_upmeso_Coelodendridae |> summary()

AIC(full_upmeso_Coelodendridae, red_upmeso_Coelodendridae)
BIC(full_upmeso_Coelodendridae, red_upmeso_Coelodendridae)

# |- LoMeso -----------------------------

full_lomeso_Coelodendridae <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(Si, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Coelodendridae, method = 'ML', select = T 
)
full_lomeso_Coelodendridae |> summary()

red_lomeso_Coelodendridae <- gam(
  intg ~ 
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Coelodendridae, method = 'ML', select = T 
)
red_lomeso_Coelodendridae |> summary()
AIC(full_lomeso_Coelodendridae, red_lomeso_Coelodendridae)
BIC(full_lomeso_Coelodendridae, red_lomeso_Coelodendridae)


###
# Collodaria #####
###

# full model
full_epi_Collodaria <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Bact_enumb, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$epi$Collodaria, method = 'ML', select = T
)
full_epi_Collodaria |>  summary()

# reduced model
red_epi_Collodaria <- gam(
  intg ~  s(sal, k = 6) + s(o2, k = 6) + 
    s(Bact_enumb, k = 6) ,
  data = taxa_rhiz$epi$Collodaria, method = 'ML', select = T
)
red_epi_Collodaria |>  summary()

AIC(full_epi_Collodaria, red_epi_Collodaria)
BIC(full_epi_Collodaria, red_epi_Collodaria)

###
# Foraminifera #############
###

# |- Epipelagic -----------------------

# full model 
full_epi_Foraminifera <- gam(
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Bact_enumb, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6),
  data = taxa_rhiz$epi$Foraminifera, method = 'ML', select = T
)
full_epi_Foraminifera |>  summary()

# reduced model
red_epi_Foraminifera <- gam(
  intg ~ s(temp, k = 6) + s(o2, k = 6) + 
    s(avg_fbn_200, k = 6) +
    s(par_conc, k = 6),
  data = taxa_rhiz$epi$Foraminifera, method = 'ML', select = T
)
red_epi_Foraminifera |>  summary()
AIC(red_epi_Foraminifera, full_epi_Foraminifera)
BIC(red_epi_Foraminifera, full_epi_Foraminifera)


###
# Plots #####
###


# |- Plot extracting function --------------------------------

smooth_predictor <- function(model, data) {
  
  # set up link function
  ilink = family(model)$linkinv
  
  # get variable predictor names
  vars <- names(model$var.summary)
  
  ## need to replace newdata with plot output  and 
  # finish merging fit, se.fit, and the newdata predicited
  # don't need to save the fitting dataframes!
  var_newdata <- list()
  
  for (var in vars) {
    
    other_vars <- vars[which(vars != var)]
    
    var_newdata[[var]] <- other_vars |> 
      lapply(get_other_median, data)
    
    names(var_newdata[[var]]) <- other_vars
    
    var_newdata[[var]][[var]] <- seq(min(data[[var]], na.rm = T), max(data[[var]], na.rm = T), length.out = 1000)
    
    var_newdata[[var]] <- as.data.frame(var_newdata[[var]])
    
    pred_fit <- predict.gam(model, var_newdata[[var]], type = 'link', se.fit = T)
    
    
  }
  
  
  
  
}

# function for inside smooth predictor
get_other_median <- function(var, data) {
  median(data[[var]], na.rm = T)
}


newdata <- data.frame(
  Si = median(tot_rhiz$upmeso$Si, na.rm = T),
  avg_mass_flux_200 = median(tot_rhiz$upmeso$avg_mass_flux_200, na.rm = T),
  par_conc = seq(min(tot_rhiz$upmeso$par_conc, na.rm = T), max(tot_rhiz$upmeso$par_conc, na.rm = T), length.out = 1000)
)


ilink = family(red_upmeso_tot)$linkinv
pred <- predict.gam(red_upmeso_tot, newdata, type = 'link', se.fit = TRUE)

ggplot() +
  geom_line(aes(x = newdata$par_conc, y = ilink(pred$fit))) +
  geom_point(aes(x = tot_rhiz$upmeso$par_conc, y = tot_rhiz$upmeso$intg))+
  theme_bw()

windows()
plot(red_upmeso_tot)
