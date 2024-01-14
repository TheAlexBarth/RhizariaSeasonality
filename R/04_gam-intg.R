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

tot_rhiz <- readRDS('./data/03_all-integrated.rds')
tot_rhiz$epi <- tot_rhiz$epi[-which(tot_rhiz$epi$sal < 33),] # drop epi malfunction
taxa_rhiz <- readRDS('./data/03_taxa-integrated.rds')

# drop epi sal malfunction
for(taxa in names(taxa_rhiz$epi)) {
  taxa_rhiz$epi[[taxa]] <- taxa_rhiz$epi[[taxa]][-which(taxa_rhiz$epi[[taxa]]$sal < 33),]
}

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

# |- Save Data for plots --------------------------

# |-|- Total Data --------------------------------------

epi_total <- smooth_predictor(red_epi_total, tot_rhiz$epi)
upmeso_total <- smooth_predictor(red_upmeso_tot, tot_rhiz$upmeso)
lomeso_total <- smooth_predictor(red_lomeso_tot, tot_rhiz$lomeso)


saveRDS(list(
  epi = epi_total,
  upmeso = upmeso_total,
  lomeso = lomeso_total
), './data/04_total-gam.rds')


# |-|- Taxa Data ----------------------------------------

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
