##
# Utilities ##
##
# Function to perform backward stepwise selection on GAM model
backward_stepwise_gam <- function(formula, data, 
                                  alpha = 0.1,
                                  method = "ML", select = TRUE) {
  # Initial model
  model <- gam(formula, data = data, method = method, select = select)
  summary_model <- summary(model)
  old_rsq <- summary_model$r.sq
  
  counter = 0
  # Loop until no non-significant terms remain
  while (TRUE) {
    # Get p-values of the terms
    pvalues <- summary_model$s.table[, "p-value"]
    
    # Check if any terms are non-significant
    if (all(pvalues < alpha)) {
      break
    }
    
    # Update the formula by removing the least significant term
    terms <- attr(terms(formula), "term.labels")
    terms <- terms[-which.max(pvalues)]
    formula <- as.formula(paste("intg ~", paste(terms, collapse = " + ")))
    
    # Refit the model with the updated formula
    new_mod <- gam(formula, data = data, method = method, select = select)
    new_summary_mod <- summary(new_mod)
    
    # check to see if rsq got much worse
    if(new_summary_mod$r.sq < (summary_model$r.sq - 0.05)) {
      break
    } else if(length(terms) < 1){
      print('No significant terms')
      break
    } else {
      model = new_mod
      summary_model = new_summary_mod
      counter = counter + 1
      cat('\n')
      cat('_________________________________________')
      cat('\n')
      cat(paste0('Model update number: ', counter))
      cat('\n')
      print(summary_model)
    }
  }
  
  return(model)
}




###
# Plot colors #######
###


seasonal_scale <- c(
  `winter` = "#96e0fa",
  `summer` = "#B2182B"
)


taxa_colors <- c(
  `Acantharea` = '#DDCC77',
  `Aulacanthidae` = '#CC6677',
  `Aulosphaeridae` = "#AA4499",
  `Castanellidae` = '#44AA99',
  `Coelodendridae` = '#332288',
  `Collodaria` = '#117733',
  `Foraminifera` = '#882255',
  `Medusettidae` = '#88CCEE',
  `Tuscaroridae` = 'black',
  `Phaeodaria` = '#E69F00',
  `Rhizaria` = 'darkgrey'
)