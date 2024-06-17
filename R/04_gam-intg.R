###
# Gam on integrated data ###
###

rm(list = ls())
set.seed(3)
library(mgcv)
library(EcotaxaTools)
library(dplyr)
library(tidyr)
library(mgcViz)
source('./R/utils.R')

##
# Load in Data #########
##

tot_rhiz <- readRDS('./data/03_all-integrated.rds')
tot_rhiz$epi <- tot_rhiz$epi[-which(tot_rhiz$epi$sal < 33),] # drop epi malfunction
taxa_rhiz <- readRDS('./data/03_taxa-integrated.rds')

# drop epi sal malfunction and RFU -999
for(taxa in names(taxa_rhiz$epi)) {
  taxa_rhiz$epi[[taxa]] <- taxa_rhiz$epi[[taxa]][-which(taxa_rhiz$epi[[taxa]]$sal < 33),]
  taxa_rhiz$epi[[taxa]] <- taxa_rhiz$epi[[taxa]][-which(taxa_rhiz$epi[[taxa]]$RFU == -999),]
}

###
# Total Rhizaria ##########
###

## |- Epi pelagic -----------------------

# full model 

full_epi_form <-   intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
  s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(pp, k = 6) + s(par_conc, k = 6)


full_epi_total <- gam(full_epi_form,
  data = tot_rhiz$epi, method = 'ML', select = T
)
full_epi_total |> summary()
# reduced model
red_epi_total <- backward_stepwise_gam(full_epi_form, tot_rhiz$epi,
                                       k = 6, method = 'ML', select = T,
                                       conc_threshold = 0.8)

## |- UpMeso ----------------------------------

full_upmeso_form <-   intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
  s(Si, k = 6) +  s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)

full_upmeso_tot <- gam(full_upmeso_form,
  data = tot_rhiz$upmeso, method = 'ML', select = T 
)
full_upmeso_tot |> summary()


red_upmeso_tot <- backward_stepwise_gam(full_upmeso_form, tot_rhiz$upmeso)

## |- LoMeso -----------------------------

full_lomeso_form <-   intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
  s(Si, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)

full_lomeso_tot <- gam(full_lomeso_form,
  data = tot_rhiz$lomeso, method = 'ML', select = T 
)
full_lomeso_tot |> summary()


red_lomeso_tot <- backward_stepwise_gam(full_lomeso_form, tot_rhiz$lomeso)


##
# Acantharea #########
##

## |- Epipelagic -----------------------

# full model \

full_epi_acantharea_form <-
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
    s(Bact_enumb, k = 6) + s(NO3, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6)

full_epi_acantharea <- gam(full_epi_acantharea_form,
  data = taxa_rhiz$epi$Acantharea, method = 'ML', select = T
)
full_epi_acantharea |>  summary()

# reduced model
red_epi_acantharea <- backward_stepwise_gam(full_epi_acantharea_form, 
                                            taxa_rhiz$epi$Acantharea)

## |- UpMeso ----------------------------------

full_upmeso_acantharea_form <- 
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6)

full_upmeso_acantharea <- gam(full_upmeso_acantharea_form,
  data = taxa_rhiz$upmeso$Acantharea, method = 'ML', select = T 
)
full_upmeso_acantharea |> summary()


red_upmeso_acantharea <- backward_stepwise_gam(full_upmeso_acantharea_form, 
                                               taxa_rhiz$upmeso$Acantharea)

## |- LoMeso -----------------------------

full_lomeso_acantharea_form <- 
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
    s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
    s(pp, k = 6) + s(par_conc, k = 6)

full_lomeso_acantharea <- gam(full_lomeso_acantharea_form,
                              data = taxa_rhiz$lomeso$Acantharea, 
                              method = 'ML', select = T)
full_lomeso_acantharea |> summary()

red_lomeso_acantharea <- backward_stepwise_gam(full_lomeso_acantharea_form, 
                                               taxa_rhiz$lomeso$Acantharea)


##
# Aulacanthidae #########
##

## |- Epipelagic -----------------------

# full model 
full_epi_Aulacanthidae_form <- 
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
  s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(pp, k = 6) + s(par_conc, k = 6) 

full_epi_Aulacanthidae <- gam(full_epi_Aulacanthidae_form,
  data = taxa_rhiz$epi$Aulacanthidae, method = 'ML', select = T
)
full_epi_Aulacanthidae |>  summary()

# reduced model
red_epi_Aulacanthidae <- backward_stepwise_gam(full_epi_Aulacanthidae_form, 
                                               taxa_rhiz$epi$Aulacanthidae)


## |- UpMeso ----------------------------------

full_upmeso_Aulacanthidae_form <- 
  intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
  s(Si, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)


full_upmeso_Aulacanthidae <- gam(full_upmeso_Aulacanthidae_form,
  data = taxa_rhiz$upmeso$Aulacanthidae, method = 'ML', select = T 
)
full_upmeso_Aulacanthidae |> summary()

red_upmeso_Aulacanthidae <- backward_stepwise_gam(full_upmeso_Aulacanthidae_form, 
                                                  taxa_rhiz$upmeso$Aulacanthidae)


## |- LoMeso -----------------------------

full_lomeso_Aulacanthidae_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
  s(Si, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)


full_lomeso_Aulacanthidae <- gam(full_lomeso_Aulacanthidae_form,
  data = taxa_rhiz$lomeso$Aulacanthidae, method = 'ML', select = T 
)
full_lomeso_Aulacanthidae |> summary()

red_lomeso_Aulacanthidae <- backward_stepwise_gam(full_lomeso_Aulacanthidae_form, 
                                                  taxa_rhiz$lomeso$Aulacanthidae)

##
# Aulosphaeridae #########
##

## |- Epipelagic -----------------------

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


## |- UpMeso ----------------------------------

full_upmeso_Aulosphaeridae_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) +
  s(Si, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)

full_upmeso_Aulosphaeridae <- gam(full_upmeso_Aulosphaeridae_form,
  data = taxa_rhiz$upmeso$Aulosphaeridae, method = 'ML', select = T 
)
full_upmeso_Aulosphaeridae |> summary()

# 
# red_upmeso_Aulosphaeridae <- backward_stepwise_gam(full_upmeso_Aulosphaeridae_form, 
#                                                   taxa_rhiz$upmeso$Aulosphaeridae)
###! backwards stepwise throws error due to multiple removals just to parconc
red_upmeso_Aulosphaeridae <- gam(
  intg ~ 
    s(par_conc, k = 6),
  data = taxa_rhiz$upmeso$Aulosphaeridae, method = 'ML', select = T 
)

## |- LoMeso -----------------------------
full_lomeso_Aulosphaeridae_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + 
  s(o2, k = 6) + s(Si, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)


full_lomeso_Aulosphaeridae <- gam(full_lomeso_Aulosphaeridae_form,
  data = taxa_rhiz$lomeso$Aulosphaeridae, method = 'ML', select = T 
)
full_lomeso_Aulosphaeridae |> summary()

# red_lomeso_Aulosphaeridae <- backward_stepwise_gam(full_lomeso_Aulosphaeridae_form, 
#                                                   taxa_rhiz$lomeso$Aulosphaeridae)

# also has error due to only par conc:
red_lomeso_Aulosphaeridae <- gam(
  intg ~ 
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Aulosphaeridae, method = 'ML', select = T 
)

###
# Castanellidae #####
###

# full model
full_epi_Castanellidae_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + 
  s(RFU, k = 6) + s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(pp, k = 6) + s(par_conc, k = 6)


full_epi_Castanellidae <- gam(full_epi_Castanellidae_form,
                              data = taxa_rhiz$epi$Castanellidae, 
                              method = 'ML', select = T)
full_epi_Castanellidae |>  summary()

# reduced model
red_epi_Castanellidae <- backward_stepwise_gam(full_epi_Castanellidae_form, 
                                               taxa_rhiz$epi$Castanellidae)

###
# Coelodendridae ###########
###


## |- UpMeso ----------------------------------

full_upmeso_Coelodendridae_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + 
  s(o2, k = 6) + s(Si, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)

full_upmeso_Coelodendridae <- gam(full_upmeso_Coelodendridae_form,
  data = taxa_rhiz$upmeso$Coelodendridae, method = 'ML', select = T)

full_upmeso_Coelodendridae |> summary()


red_upmeso_Coelodendridae <- backward_stepwise_gam(full_upmeso_Coelodendridae_form, 
                                                   taxa_rhiz$upmeso$Coelodendridae)

## |- LoMeso -----------------------------

full_lomeso_Coelodendridae_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + 
  s(o2, k = 6) + s(Si, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(par_conc, k = 6)

full_lomeso_Coelodendridae <- gam(full_lomeso_Coelodendridae_form,
  data = taxa_rhiz$lomeso$Coelodendridae, method = 'ML', select = T 
)
full_lomeso_Coelodendridae |> summary()

# red_lomeso_Coelodendridae <- backward_stepwise_gam(full_lomeso_Coelodendridae_form, 
#                                                    taxa_rhiz$lomeso$Coelodendridae)
# 
# errors due to one term:
red_lomeso_Coelodendridae <- gam(
  intg ~ 
    s(par_conc, k = 6),
  data = taxa_rhiz$lomeso$Coelodendridae, method = 'ML', select = T 
)

###
# Collodaria #####
###

# full model
full_epi_Collodaria_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + 
  s(RFU, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(pp, k = 6) + s(par_conc, k = 6)

full_epi_Collodaria <- gam(full_epi_Collodaria_form,
                           data = taxa_rhiz$epi$Collodaria, 
                           method = 'ML', select = T)
full_epi_Collodaria |>  summary()

# reduced model
red_epi_Collodaria <- backward_stepwise_gam(full_epi_Collodaria_form, 
                                            taxa_rhiz$epi$Collodaria)

###
# Foraminifera #############
###

## |- Epipelagic -----------------------

# full model 
full_epi_Formainifera_form <- intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + 
  s(RFU, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
  s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
  s(pp, k = 6) + s(par_conc, k = 6)


full_epi_Foraminifera <- gam(full_epi_Formainifera_form,
  data = taxa_rhiz$epi$Foraminifera, method = 'ML', select = T
)
full_epi_Foraminifera |>  summary()

# reduced model
red_epi_Foraminifera <- backward_stepwise_gam(full_epi_Formainifera_form, 
                                              taxa_rhiz$epi$Foraminifera)



####
# Model Summaries #############
####


red_epi_total |>  summary()
red_upmeso_tot |> summary()
red_lomeso_tot |> summary()
red_epi_acantharea |>  summary()
red_upmeso_acantharea |> summary()
red_lomeso_acantharea |> summary()
red_epi_Aulacanthidae |>  summary()
red_upmeso_Aulacanthidae |> summary()
red_lomeso_Aulacanthidae |> summary()
red_upmeso_Aulosphaeridae |> summary()
red_lomeso_Aulosphaeridae |> summary()
red_epi_Castanellidae |>  summary()
red_upmeso_Coelodendridae |> summary()
red_lomeso_Coelodendridae |> summary()
red_epi_Collodaria |>  summary()
red_epi_Foraminifera |>  summary()


## |- Save Models --------------------------
saveRDS(
  list(
      tot_epi = red_epi_total,
      tot_upmeso = red_upmeso_tot,
      tot_lomeso = red_lomeso_tot,
      acantharea_epi = red_epi_acantharea,
      acantharea_upmeso = red_upmeso_acantharea,
      acantharea_lomeso = red_lomeso_acantharea,
      aulacanthidae_epi = red_epi_Aulacanthidae,
      aulacanthidae_upmeso = red_upmeso_Aulacanthidae,
      aulacanthidae_lomeso = red_lomeso_Aulacanthidae,
      aulosphaeridae_upmeso = red_upmeso_Aulosphaeridae,
      aulosphaeridae_lomeso = red_lomeso_Aulosphaeridae,
      castanellidae_epi = red_epi_Castanellidae,
      coelodendridae_upmeso = red_upmeso_Coelodendridae,
      coelodendridae_lomeso = red_lomeso_Coelodendridae,
      collodaria_epi = red_epi_Collodaria,
      foraminifera_epi = red_epi_Foraminifera
    ),
  './data/04_gam-intg.rds'
)





###
# Plots #####
###


## |- Plot extracting function --------------------------------

# function for inside smooth predictor
get_other_median <- function(var, data) {
  median(data[[var]], na.rm = T)
}

smooth_predictor <- function(model, data) {
  
  # set up link function
  ilink = family(model)$linkinv
  
  # get variable predictor names
  vars <- names(model$var.summary)
  
  ## need to replace newdata with plot output  and 
  # finish merging fit, se.fit, and the newdata predicited
  # don't need to save the fitting dataframes!
  
  plot_output <- list()
  
  for (var in vars) {
    
    other_vars <- vars[which(vars != var)]
    
    var_newdata<- other_vars |> 
      lapply(get_other_median, data)
    
    names(var_newdata) <- other_vars
    
    var_newdata[[var]] <- seq(min(data[[var]], na.rm = T), max(data[[var]], na.rm = T), length.out = 1000)
    
    var_newdata <- as.data.frame(var_newdata)
    
    pred_fit <- predict.gam(model, var_newdata, type = 'link', se.fit = T)
    
    plot_output[[var]] <- list()
    plot_output[[var]]$pred <- data.frame(
      fit = ilink(pred_fit$fit),
      low = ilink(pred_fit$fit - (1.96*pred_fit$se.fit)),
      high = ilink(pred_fit$fit + (1.96 * pred_fit$se.fit))
    )
    plot_output[[var]]$pred[[var]] <- var_newdata[[var]]
    
    plot_output[[var]]$obs <- data.frame(
      intg = data$intg,
      month = data$month
    )
    plot_output[[var]]$obs[[var]] <- data[[var]]
    
  }
  
  return(plot_output)
}

## |- Save Data for plots --------------------------

## |-|- Total Data --------------------------------------

epi_total <- smooth_predictor(red_epi_total, tot_rhiz$epi)
upmeso_total <- smooth_predictor(red_upmeso_tot, tot_rhiz$upmeso)
lomeso_total <- smooth_predictor(red_lomeso_tot, tot_rhiz$lomeso)


saveRDS(list(
  epi = epi_total,
  upmeso = upmeso_total,
  lomeso = lomeso_total
), './data/04_total-gam.rds')


## |-|- Taxa Data ----------------------------------------

acantharea_epi <- smooth_predictor(red_epi_acantharea, taxa_rhiz$epi$Acantharea)
acantharea_upmeso <- smooth_predictor(red_upmeso_acantharea, taxa_rhiz$upmeso$Acantharea)
acantharea_lomeso <- smooth_predictor(red_lomeso_acantharea, taxa_rhiz$lomeso$Acantharea)
saveRDS(
  list(
    epi = acantharea_epi,
    upmeso = acantharea_upmeso,
    lomeso = acantharea_lomeso
  ),
  './data/04_acanth-gam.rds'
)

aulacanthidae_epi <- smooth_predictor(red_epi_Aulacanthidae, taxa_rhiz$epi$Aulacanthidae)
aulacanthidae_upmeso <- smooth_predictor(red_upmeso_Aulacanthidae, taxa_rhiz$upmeso$Aulacanthidae)
aulacanthidae_lomeso <- smooth_predictor(red_lomeso_Aulacanthidae, taxa_rhiz$lomeso$Aulacanthidae)
saveRDS(
  list(
    epi = aulacanthidae_epi,
    upmeso = aulacanthidae_upmeso,
    lomeso = aulacanthidae_lomeso
  ),
  './data/04_aulacanthidae-gam.rds'
)

aulosphaeridae_upmeso <- smooth_predictor(red_upmeso_Aulosphaeridae, taxa_rhiz$upmeso$Aulosphaeridae)
aulosphaeridae_lomeso <- smooth_predictor(red_lomeso_Aulosphaeridae, taxa_rhiz$lomeso$Aulosphaeridae)
saveRDS(
  list(
    upmeso = aulosphaeridae_upmeso,
    lomeso = aulosphaeridae_lomeso
  ),
  './data/04_aulosphaeridae-gam.rds'
)

castanellidae_epi <- smooth_predictor(red_epi_Castanellidae, taxa_rhiz$epi$Castanellidae)
saveRDS(
  list(
   epi = castanellidae_epi
  ),
  './data/04_castanellidae-gam.rds'
)

coelodendridae_upmeso <- smooth_predictor(red_upmeso_Coelodendridae, taxa_rhiz$upmeso$Coelodendridae)
coelodendridae_lomeso <- smooth_predictor(red_lomeso_Coelodendridae, taxa_rhiz$lomeso$Coelodendridae)
saveRDS(
  list(
    upmeso = coelodendridae_upmeso,
    lomeso = coelodendridae_lomeso
  ),
  './data/04_coelodendridae-gam.rds'
)


collodaria_epi <- smooth_predictor(red_epi_Collodaria, taxa_rhiz$epi$Collodaria)
saveRDS(
  list(
    epi = collodaria_epi
  ),
  './data/04_collodaria-gam.rds'
)


foraminifera_epi <- smooth_predictor(red_epi_Foraminifera, taxa_rhiz$epi$Foraminifera)
saveRDS(
  list(
    epi = foraminifera_epi
  ),
  './data/04_foram-gam.rds'
)
