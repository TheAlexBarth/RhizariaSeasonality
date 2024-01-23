###
# Plot percent Contribution ######
###

rm(list = ls())

library(EcotaxaTools)
library(dplyr)
library(ggplot2)

uvp_meta <- readRDS('./data/00_zoop-uvp.rds')$meta


all_den <- readRDS("./data/s02_all-taxa-density.rds")

all_df <- list_to_tib(all_den, 'profileid')
all_df$zone <- NA
for(r in 1:nrow(all_df)) {
  if(all_df$max_d[r] <= 200) {
    all_df$zone[r] <- 'epi'
  } else if(all_df$max_d[r] <= 500) {
    all_df$zone[r] <- 'upmeso'
  } else {
    all_df$zone[r] <- 'lomeso'
  }
}


all_avg <- all_df |> 
  left_join(
    uvp_meta |> 
      select(profileid, cruise_id, sampledate)
  ) |> 
  group_by(group, zone, cruise_id) |> 
  summarize(avg = mean(conc_m3),
            sd = sd(conc_m3))


ggplot(all_avg) +
  geom_bar(aes(x = cruise_id, y = avg, fill = group),
           stat = 'identity', position = 'dodge')+
  geom_errorbar(aes(x = cruise_id,
                    ymin = avg,
                    ymax = avg + sd,
                    group = group),
                stat = 'identity', position = 'dodge') +
  facet_wrap(.~zone, ncol = 1)+
  theme_bw()

ggsave('./output/Supplement/s02-total-contrib.pdf',
       width = 190, dpi = 600, units = 'mm')
