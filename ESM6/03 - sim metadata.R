# # --------------------------------------------
# # -- INFO needed to calc stats
# # --------------------------------------------
SIM_DIR_NAME <- 'C:/Users/esokol/Box/PROJECTS/B235/Metacomm Sims/R_sims_v2_10000_reps/SIM_RESULTS'

# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
options(stringsAsFactors = FALSE)

require(tidyverse)

# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
# SIM_DIR_NAME<-paste('SIM',scenario.ID,sep='_')  # Directory name

sim.file.list<-list.files(SIM_DIR_NAME)

RDATA.list_ALL <- sim.file.list[grep('.rda',sim.file.list)]

set.seed(1234)
# RDATA.list <- RDATA.list_ALL %>% sample(1000) #subsample files
RDATA.list <- RDATA.list_ALL



# -- loop to do calcs for each sim

d_results <- tibble()
for (i in 1:length(RDATA.list)){
  i.RDATA <- RDATA.list[i]
  d_results_tmp <- tibble()
  try({
    i.RDATA.filename<-paste(SIM_DIR_NAME,i.RDATA,sep='/')
    
    load(i.RDATA.filename)
    
    sim.ID <- sim.result$sim.result.name
    scenario.ID <- sim.result$scenario.ID
    
    d_site_info <- sim.result$landscape$site.info
    
    d_results_tmp <- tibble(
      scenario_id = scenario.ID,
      sim_file_name = i.RDATA.filename,
      sim_id = sim.ID,
      m_black_mats = d_site_info$m_values[d_site_info$mat_type == 'black'][1],
      m_other = d_site_info$m_values[d_site_info$mat_type != 'black'][1],
      proportion_black_mats = sum(d_site_info$mat_type == 'black')/length(d_site_info$mat_type),
      proportion_orange_mats = sum(d_site_info$mat_type == 'orange')/length(d_site_info$mat_type),
      proportion_red_mats = sum(d_site_info$mat_type == 'red')/length(d_site_info$mat_type),
      proportion_green_mats = sum(d_site_info$mat_type == 'green')/length(d_site_info$mat_type),
      proportion_biomass_black_mats = sum(as.numeric(d_site_info$mat_type == 'black') * d_site_info$JL) / sum(d_site_info$JL),
      proportion_biomass_orange_mats = sum(as.numeric(d_site_info$mat_type == 'orange') * d_site_info$JL) / sum(d_site_info$JL),
      proportion_biomass_red_mats = sum(as.numeric(d_site_info$mat_type == 'red') * d_site_info$JL) / sum(d_site_info$JL),
      proportion_biomass_green_mats = sum(as.numeric(d_site_info$mat_type == 'green') * d_site_info$JL) / sum(d_site_info$JL)
    )
    
  })
  
  d_results <- bind_rows(d_results, d_results_tmp)
  
  print(paste(i, 'out of', length(RDATA.list)))
}

d_results <- d_results %>%
  separate(sim_id, into = c('sim','metacomm','niche_scaling','m','nu','w','rep','date','time','date_time'),
           sep = '_',
           remove = FALSE,
           extra = 'merge') %>% 
  mutate(niche_scaling = gsub('nicheScaling\\-','',niche_scaling),
         nu = gsub('nu\\-','',nu),
         w = gsub('w\\-','',w),
         m_ratio_black_over_other = m_black_mats / m_other) %>%
  select(-m, -date, -time, -date_time, -rep, -sim, - metacomm)

# write out useful metadata
readr::write_csv(
  d_results, 
  path = paste('d_sim_metadata_',scenario.ID,'.csv',sep=''))
