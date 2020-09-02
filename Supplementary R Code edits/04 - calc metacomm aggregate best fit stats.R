#' R Script 04 for:
#' Sokol ER, Barrett JE, Kohler TJ, McKnight DM, Salvatore MR and Stanish LF (2020)
#'  Evaluating Alternative Metacommunity Hypotheses for Diatoms in the McMurdo Dry
#'  Valleys Using Simulations and Remote Sensing Data. Front. Ecol. Evol. 8:521668. 
#'  doi: 10.3389/fevo.2020.521668
#'
#' @author Eric R. Sokol \email{esokol@battelleecology.org}
#' 
#' This R script calculates alpha, beta, and gamma diversity for each simulation outcome
#' at timesteps 1 and 100. Then calculates Chi-squared fit statistics as described in 
#' the article cited above.



# clear out workspace
rm(list=ls())
gc()


# set path to directory with sim outputs
# SIM_DIR_NAME <- MY_PATH_TO_SIMULATION_OUTPUT_DIRECTORY


# set options
options(stringsAsFactors = FALSE)


# load required packaqes
require(vegetarian)
require(tidyverse)


# set paths for reading and writing data
# you may need to alter your path to the directory with this table
read_file_path <- "Supplementary R Code"
write_file_path <- "Supplementary R Code"


# read in metadata file
d_sim_metadata <- read_csv(paste0(read_file_path,"/","d_sim_metadata_streamMetacommunities.csv"))


# modify sim_id to create sim filenames that can be read into R env
RDATA.list <- d_sim_metadata$sim_id %>% paste0(".rda")



##################################
##################################
##################################
##################################

# loop to calculate biodiversity summary statistics for first and final time step
# of each simulation.Summary stats include alpha, beta, and gamma diversity



# i.RDATA <- RDATA.list[1]


# wrapper function that can be called for parallel processing
fnx <- function(i.RDATA,
                time_steps_to_use = c(1,100),
                SIM_DIR_NAME = "SIM_RESULTS"){

  d_results <- data.frame()
  
  try({
    sim.result <- NULL
    i.RDATA.filename<-paste(SIM_DIR_NAME,i.RDATA,sep="/")
    
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
                    return(vegetarian::d(abundances = x, lev = "alpha", q = 1))
                  }) %>% unlist(),
      beta_q1 = map(data_wide,
                     function(x){
                       x <- x %>% as.data.frame()
                       return(vegetarian::d(abundances = x, lev = "beta", q = 1))
                     }) %>% unlist(),
      gamma_q1 = map(data_wide,
                    function(x){
                      x <- x %>% as.data.frame()
                      return(vegetarian::d(abundances = x, lev = "gamma", q = 1))
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
      note = "error in div calc")
    print(paste0("ERROR in div calc for ",i.RDATA))
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


# write out results
readr::write_csv(d_metacomm_stats, path = paste0(write_file_path,"/d_metacomm_aggregate_stats_",Sys.Date(), ".csv"))


# calc fit stats
sd_alpha = sd(d_metacomm_stats$alpha_q1)
sd_beta = sd(d_metacomm_stats$beta_q1)
sd_gamma = sd(d_metacomm_stats$gamma_q1)


# make data frame with summary stats
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
write_csv(d_metacomm_stats, paste0(write_file_path,"/d_metacomm_fit_stats_aggregate_" ,Sys.Date(), ".csv"))


