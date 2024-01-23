###
# Calc Densities for Rhizaria and other groups ####
###


uvp_data <- readRDS("./data/00_zoop-uvp.rds")

####
# Re-naming zoops ########
####

taxa_names <- c('Rhizaria', 'Copepoda', 'living')

uvp_data <- uvp_data |> 
  add_zoo(func = names_to, col_name = 'name',
          new_names = taxa_names)

###
# Calculate Binned Abundances #############
###

# |- Rhizaria Only Calcs ----------------------------------------
all_density <- uvp_data |> 
  uvp_zoo_conc(breaks = seq(0,1000,25))

# |- Add 0 observations for certain UVP taxa --------------------
all_possible_taxa <- unique(unlist(lapply(all_density, function(x) x$group)))

# Function to merge each dataframe in the list with the template dataframe

fill_missing_species <- function(df) {
  if(all(all_possible_taxa %in% unique(df$group))) {
    return(df)
  } else {
    merge_df <- expand.grid(
      db = unique(df$db),
      group = all_possible_taxa
    )
    
    df <- df |> 
      merge(merge_df, by = c('group', 'db'), all.y = T)
    
    df[is.na(df)] <- 0
    return(df)
  }
}

corrected_densities <- all_density |>
  lapply(fill_missing_species) |>
  lapply(bin_format)


saveRDS(corrected_densities,"./data/s02_all-taxa-density.rds")
