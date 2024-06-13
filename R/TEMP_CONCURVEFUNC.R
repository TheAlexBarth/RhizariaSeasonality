library(mgcv)

# Function to fit GAM and check for concurvity
fit_gam_check_concurvity <- function(formula, data, k = 6, method = "ML", 
                                     select = TRUE, threshold = 0.8) {
  while (TRUE) {
    # Fit the model
    model <- gam(formula, data = data, method = method, select = select)
    
    # Check concurvity
    conc <- concurvity(model, full = FALSE)$estimate
    #format the matrix
    conc <- conc[2:nrow(conc),]
    conc <- conc[,2:ncol(conc)]
    conc[lower.tri(conc, diag = T)] <- NA
    lower.tri(conc)
    
    # Find pairs with concurvity above the threshold
    # If no high concurvity pairs found, break the loop
    if (all(conc < threshold)) {
      break
    } else {
      max_conc <- max(conc, na.rm = T)
      max_idx <- which(conc == max_conc, arr.ind = T)
      
      if(nrow(max_idx) > 1) {
        cat(paste0('\n', 'There are more than one max concurvity scores',
                   '\n', 'They are for terms: ', rownames(max_idx), '\n',
                   'I didn\'t plan for this, fix funciton in utils or self evaluate'))
      } else {
        term1 <- rownames(conc)[max_idx[1]] |> add_ks(k)
        term2 <- colnames(conc)[max_idx[2]]
        
        ### STOP HERE TO FIX issue with adding ks to the model terms
        
        
        
      }
    }
    # Identify terms in the pair with the highest concurvity
    term1 <- rownames(high_conc_pairs)[1]
    term2 <- colnames(high_conc_pairs)[1]
    
    # Create formulas with each term removed
    terms <- attr(terms(formula), "term.labels")
    formula1 <- as.formula(paste("intg ~", paste(setdiff(terms, term1), collapse = " + ")))
    formula2 <- as.formula(paste("intg ~", paste(setdiff(terms, term2), collapse = " + ")))
    
    # Fit models with each term removed
    model1 <- gam(formula1, data = data, method = method, select = select)
    model2 <- gam(formula2, data = data, method = method, select = select)
    
    # Compare models (e.g., using AIC)
    aic1 <- AIC(model1)
    aic2 <- AIC(model2)
    
    # Select the model with the lower AIC
    if (aic1 < aic2) {
      formula <- formula1
    } else {
      formula <- formula2
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


# # Initial formula
# initial_formula <- intg ~ s(temp, k = 6) + s(sal, k = 6) + s(o2, k = 6) + s(RFU, k = 6) +
#   s(Si, k = 6) + s(Bact_enumb, k = 6) + s(NO3, k = 6) +
#   s(avg_mass_flux_200, k = 6) + s(avg_fbc_200, k = 6) + s(avg_fbn_200, k = 6) +
#   s(pp, k = 6) + s(small, k = 6) + s(large, k = 6)
# 
# # Data
# data <- tot_rhiz$epi
# 
# # Run the function to fit the model and check concurvity
# final_model <- fit_gam_check_concurvity(initial_formula, data)
# 
# # Print the summary of the final model
# summary(final_model)
