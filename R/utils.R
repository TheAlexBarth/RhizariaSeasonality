##
# Utilities ##
##



#' Function to perform backward stepwise selection on GAM model
#' 
#' This function performs backward stepwise selection on a GAM model by
#' iteratively removing the least significant term until all terms are
#' significant. The function also checks for concurvity in the model and
#' refits the model if concurvity is detected.
#' 
#' @param formula a formula object specifying the model
#' @param data a data frame containing the data
#' @param alpha the significance level for removing terms
#' @param k the number of knots to use in the model
#' @param method the method to use for fitting the model
#' @param select whether to use the select argument in the gam function
#' @param conc_threshold the threshold for concurvity
backward_stepwise_gam <- function(formula, data, 
                                  alpha = 0.1, k = 6,
                                  method = "ML", select = TRUE,
                                  conc_threshold = 0.8) {
  
  # Initial model passed to concurvity check:
  model <- fit_gam_check_concurvity(formula, data, 
                                    k = k, 
                                    method = method, 
                                    select = select,
                                    threshold = conc_threshold)
  
  summary_model <- summary(model)
  old_rsq <- summary_model$r.sq
  
  counter = 0
  # Loop until no non-significant terms remain
  while (TRUE) {
    # Get p-values of the terms
    pvalues <- summary_model$s.table[, "p-value"]
    
    # Check if any terms are non-significant
    if (all(pvalues < alpha)) {
      cat(paste0('\n', 'Breaking due to alpha threshold'))
      break
    }
    
    # Update the formula by removing the least significant term
    formula <- formula(model)
    terms <- attr(terms(formula), "term.labels")
    terms <- terms[-which.max(pvalues)]
    formula <- as.formula(paste("intg ~", paste(terms, collapse = " + ")))
    
    # Refit the model with the updated formula
    # this is the proposal model
    # need to concurvity check the new model!
    new_mod <- fit_gam_check_concurvity(formula, data, 
                                        k = k, 
                                        method = method, 
                                        select = select)
    new_summary_mod <- summary(new_mod)
    
    # check to see if rsq got much worse
    if(new_summary_mod$r.sq < (summary_model$r.sq - 0.05)) {
      cat('\n', 'breaking due to rsq drop')
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


# Function to fit GAM and check for concurvity
fit_gam_check_concurvity <- function(formula, data, k = 6, method = "ML", 
                                     select = TRUE, threshold = 0.8) {
  
  model <- gam(formula, data = data, method = method, select = select)
  while (TRUE) {
    
    # Check concurvity
    conc <- concurvity(model, full = FALSE)$estimate
    #format the matrix to remove para
    conc <- conc[2:nrow(conc),]
    conc <- conc[,2:ncol(conc)]
    conc[lower.tri(conc, diag = T)] <- NA
    
    # Find pairs with concurvity above the threshold
    # If no high concurvity pairs found, break the loop
    if (all(conc[!is.na(conc)] < threshold)) {
      cat('\n', 'Concurvity check passed!')
      break
    } else {
      max_conc <- max(conc, na.rm = T)
      max_idx <- which(conc == max_conc, arr.ind = T)
      
      if(nrow(max_idx) > 1) {
        cat(paste0('\n', 'There are more than one max concurvity scores',
                   '\n', 'They are for terms: ', rownames(max_idx), '\n',
                   'I didn\'t plan for this, fix funciton in utils or self evaluate'))
      } else {
        term1 <- rownames(conc)[max_idx[1]] |> add_ks(k = k)
        term2 <- colnames(conc)[max_idx[2]] |> add_ks(k = k)
        
        cat(paste0('\n','Concurvity found between: ', term1, ' and ', term2,
                   '\n', 'Score was: ', max_conc, '\n', 'Refitting...'))
        
        all_terms <- attr(terms(formula), 'term.labels')
        
        mod1_formula <- paste('intg ~', 
                              paste(setdiff(all_terms, term1), 
                                    collapse = ' + ')) |> 
          as.formula()
        
        
        mod2_formula <- paste('intg ~', 
                              paste(setdiff(all_terms, term2), 
                                    collapse = '+')) |> 
          as.formula()
        
        ## Check and compare models
        mod1 <- gam(mod1_formula, data = data, method = method, select = select)
        mod2 <- gam(mod2_formula, data = data, method = method, select = select)
        
        aic1 <- AIC(mod1)
        aic2 <- AIC(mod2)
        
        if(aic1 < aic2) {
          cat(paste0('\n', 'Dropping: ', term1, ' due to concurvity'))
          model <- mod1
          formula <- mod1_formula
        } else if(aic2 < aic1){
          cat(paste0('\n', 'Dropping: ', term2, ' due to concurvity'))
          model <- mod2
          formula <- mod2_formula
        } else {
          model <- sample(c(mod1, mod2), 1)
          if(model == mod1) {
            cat(paste0('\n', 'Dropping: ', term2, ' due to concurvity at random'))
            formula <- mod1_formula
          } else {
            cat(paste0('\n', 'Dropping: ', term1, ' due to concurvity at random'))
            formula <- mod2_formula
          }
        }
      }
    }
  }
  
  return(model)
}


#' Add a k value to the term functions
#'
#' There is a probably better option but this was needed to make it 
#' run smoothly with the code I wrote
#' 
#' @param term a term for a gam model, MUST BE formatted as s(term)
#' @param k the nubmer of knots
add_ks <- function(term, k) {
  sub(')', paste0(', k = ',k,')'),term)
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