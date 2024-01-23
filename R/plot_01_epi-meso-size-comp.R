###
# Are acantharea rhizaria smaller?
###



rm(list = ls())

library(ggplot2)
library(ggpubr)
library(dplyr)
library(EcotaxaTools)

###
# Set up ######################
###

uvp_data <- readRDS("./data/00_zoop-uvp.rds")

####
# Re-naming zoops ########
####

taxa_names <- c('Rhizaria', 'Copepoda',
                'Phaeodaria','Aulacanthidae','Coelodendridae',
                'Aulosphaeridae','Acantharea',
                'Castanellidae','Collodaria',
                'Medusettidae','Foraminifera',
                'Tuscaroridae','living')

uvp_data <- uvp_data |> 
  add_zoo(func = names_to, col_name = 'name',
          new_names = taxa_names)

###
# Rhiz only #########
###

rhiz_only <- uvp_data |> 
  mod_zoo(func = 'names_keep', keep_names = 'Rhizaria',
          keep_children = T)


# remove data.frames with 0's
no_rows <- NULL
for(name in names(rhiz_only$zoo_files)) {
  if(nrow(rhiz_only$zoo_files[[name]]) == 0) {
    no_rows <- c(no_rows, name)
  }
}

rhiz_only$par_files <- rhiz_only$par_files[-which(names(rhiz_only$par_files) %in% no_rows)]
rhiz_only$zoo_files <- rhiz_only$zoo_files[-which(names(rhiz_only$zoo_files) %in% no_rows)]
rhiz_only$meta <- rhiz_only$meta[-which(rhiz_only$meta$profileid %in% no_rows),]

rhiz_sizes <- rhiz_only$zoo_files |> 
  list_to_tib('profileid')
rhiz_sizes$esd <- rhiz_sizes$esd*unique(rhiz_only$meta$acq_pixel)


####
# Size by depth group #
####

rhiz_sizes$zone <- NA
for(r in 1:nrow(rhiz_sizes)) {
  if(rhiz_sizes$depth_including_offset[r] <= 200) {
    rhiz_sizes$zone[r] <- 'epi'
  } else{
    rhiz_sizes$zone[r] <- 'meso'
  }
}

# |- filter to more than 100 obs --------------

rhiz_sizes <- rhiz_sizes |> 
  filter(name %in% c('Acantharea', 'Aulacanthidae', 'Foraminifera', 'Aulosphaeridae'))

# |- Plot sizes -------------------------

ggplot(rhiz_sizes,aes(x = name, y = esd, fill = zone, color = zone)) +
  geom_point(position = position_jitterdodge(0.2), alpha = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  scale_fill_manual(values = c(
    `epi` = '#5acadb',
    `meso` = '#014d59'
  )) +
  scale_color_manual(values = c(
    `epi` = '#5acadb',
    `meso` = '#014d59'
  )) +
  theme_pubr() +
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave('./output/09_size-comp/size-comp.pdf',width = 85, units = 'mm')


# |- Stats ----------------------
for(name in unique(rhiz_sizes$name)) {
  print(name)
  loop_data <- rhiz_sizes[which(rhiz_sizes$name == name),]
  wilcox.test(esd ~ zone, loop_data) |> 
    print()
  mean(loop_data$esd[loop_data$zone == 'epi']) |> print()
  mean(loop_data$esd[loop_data$zone == 'meso']) |> print()
  print(mean(loop_data$esd[loop_data$zone == 'epi']) - mean(loop_data$esd[loop_data$zone == 'meso']))
}
