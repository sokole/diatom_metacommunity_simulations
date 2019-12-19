rm(list=ls())
gc()

# # --------------------------------------------
# # -- INFO needed to calc stats
# # --------------------------------------------
# scenario.ID<-'DISPERSAL_EXAMPLES'
SIM_DIR_NAME <- 'C:/Users/esokol/Box/PROJECTS/B235/Metacomm Sims/R_sims_v2_10000_reps/SIM_RESULTS'

# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
options(stringsAsFactors = FALSE)

require(vegetarian)
require(tidyverse)

#######################################
#######################################
#######################################
#######################################
# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
# --------------------------------------------
# SIM_DIR_NAME<-paste('SIM',scenario.ID,sep='_')  # Directory name

d_sim_metadata <- read_csv('d_sim_metadata_streamMetacommunities.csv')
RDATA.list <- d_sim_metadata$sim_id %>% paste0('.rda')



##################################
##################################
##################################
##################################
# i.RDATA <- RDATA.list[1]

fnx <- function(i.RDATA,
                time_steps_to_use = c(1,100),
                SIM_DIR_NAME = 'SIM_RESULTS'){

  d_results <- data.frame()
  
  try({
    sim.result <- NULL
    i.RDATA.filename<-paste(SIM_DIR_NAME,i.RDATA,sep='/')
    
    load(i.RDATA.filename)
    
    sim.ID <- sim.result$sim.result.name
    scenario.ID <- sim.result$scenario.ID
    
    d_site_info <- sim.result$landscape$site.info
    
    # J.long <- sim.result$J.long
    J.long <- sim.result$J.long %>% filter(timestep %in% time_steps_to_use)
    J.long.nested <- J.long %>% nest(-timestep)  
    d_results_nested <- J.long.nested %>% mutate(
      data_wide = map(data,
                      function(x){ 
                        x %>% as.data.frame() %>% spread(spp, -site) %>% select(-site)
                      }),
      alpha_q1 = map(data_wide,
                  function(x){
                    x <- x %>% as.data.frame()
                    return(vegetarian::d(abundances = x, lev = 'alpha', q = 1))
                  }) %>% unlist(),
      beta_q1 = map(data_wide,
                     function(x){
                       x <- x %>% as.data.frame()
                       return(vegetarian::d(abundances = x, lev = 'beta', q = 1))
                     }) %>% unlist(),
      gamma_q1 = map(data_wide,
                    function(x){
                      x <- x %>% as.data.frame()
                      return(vegetarian::d(abundances = x, lev = 'gamma', q = 1))
                    }) %>% unlist())
    
    
    
    d_results <- data.frame(
      scenario_id = scenario.ID,
      sim_id = sim.ID,
      sim_filename = i.RDATA,
      d_results_nested %>% select(-data, -data_wide)) 
    
  })
  
  if (nrow(d_results) == 0){
    d_results <- data.frame(
      scenario_id = scenario.ID,
      sim_id = sim.ID,
      sim_filename = i.RDATA,
      note = 'error in div calc')
    print(paste0('ERROR in div calc for ',i.RDATA))
  }else{
    print(i.RDATA)
  }
  
  return(d_results)
  
}
##################################
##################################
##################################
##################################

# apply above fnx to each item in RDATA.list
d_metacomm_stats <- apply(
  X = RDATA.list %>% as.array(), 
  MARGIN = 1,
  FUN = fnx,
  SIM_DIR_NAME = SIM_DIR_NAME) %>% bind_rows()

# # remove duplicates
# d_metacomm_unique_obs <- d_metacomm_stats %>% 
#   select(sim_id, timestep) %>%
#   distinct()
# 
# d_metacomm_stats %>%
#   filter(sim_id == 'SIM_METACOMM_nicheScaling-4_m-0.1_nu-0.001_w-1e+05_rep-854_2019-08-29_200713_20190829_200941')
  
# write out results
readr::write_csv(d_metacomm_stats, path = paste0('d_metacomm_aggregate_stats_',Sys.Date(), '.csv'))
# d_metacomm_stats <- readr::read_csv('d_metacomm_aggregate_stats_2019-09-10.csv')


# calc fit stats
sd_alpha = sd(d_metacomm_stats$alpha_q1)
sd_beta = sd(d_metacomm_stats$beta_q1)
sd_gamma = sd(d_metacomm_stats$gamma_q1)

# d_metacomm_stats_summary <- d_metacomm_stats %>%
#   group_by(sim_id) %>%
#   summarize(n_obs = n())
# d_metacomm_stats_summary %>% filter(n_obs != 2)

d_metacomm_stats <- d_metacomm_stats %>% group_by(sim_id) %>%
  summarise_at(
    vars(alpha_q1, beta_q1, gamma_q1),
    funs(length, diff)
  ) %>% mutate(
    alpha_chi_sq = alpha_q1_diff^2 / sd_alpha,
    beta_chi_sq = beta_q1_diff^2 / sd_beta,
    gamma_chi_sq = gamma_q1_diff^2 / sd_gamma,
    chi_sq_tot = alpha_chi_sq + beta_chi_sq + gamma_chi_sq)

# write out fit stats
write_csv(d_metacomm_stats, paste0('d_metacomm_fit_stats_aggregate_' ,Sys.Date(), '.csv'))


