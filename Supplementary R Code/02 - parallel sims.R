#' R Script 02 for:
#' Sokol ER, Barrett JE, Kohler TJ, McKnight DM, Salvatore MR and Stanish LF (2020)
#'  Evaluating Alternative Metacommunity Hypotheses for Diatoms in the McMurdo Dry
#'  Valleys Using Simulations and Remote Sensing Data. Front. Ecol. Evol. 8:521668. 
#'  doi: 10.3389/fevo.2020.521668
#'
#' @author Eric R. Sokol \email{esokol@battelleecology.org}
#' 
#' This R script imports NDVI data and MCM LTER archive data of diatoms counts
#' in cyanobacteria mat cores sampled from stream mats in Fryxell Basin in Taylor
#' Valley in Antarctica, conducts an OMI analysis, and then creates 10 000 
#' metacommunity simulations using MCSim (v0.4.9). 



###############################################
# set seed should create same simulation outcomes as reported in the paper
# Note that simulations were run using R v3.6. The set.seed() function might
# produce different results in R v4+

set.seed(1234)

# select number of parallel threads to run
n_cores <- 15

# number of simulations
n_sims <- 10000

# required packages for this script
library(tidyverse)
library(parallel)
library(ade4)
library(vegan)

# Simulations used MCSim v0.4.9 with R v3.6 
# devtools::install_github('sokole/MCSim@v0.4.9')
library(MCSim)

# modify options
options(stringsAsFactors = FALSE)


#############################################################
# read in diatom community data
df_diatoms_in <- read_csv('Supplementary_Material_Data_Sheet_1_LTER_Algae_Ops_data.csv') %>%
  filter(
    !grepl('(?i)onyx',`TRANSECT NAME`) & 
      !grepl('(?i)adams',`TRANSECT NAME`) &
      !grepl('(?i)garwood', `TRANSECT NAME`) &
      !grepl('(?i)miers', `TRANSECT NAME`) &
      !is.na(`CHL-A (ug/cm2)`) &
      !is.na(`TOTAL_ANNUAL_DISCHARGE (L)`) &
      !is.na(`AFDM (mg/cm2)`)) %>%
  arrange(desc(`CHL-A (ug/cm2)`))

# only include orange and black mat records
df_diatoms_in <- df_diatoms_in %>%
  filter(grepl('(?i)orange',`MAT TYPE`) | grepl('(?i)black',`MAT TYPE`))

# group records by "region"
rows_canada_huey <- which(grepl('(?i)canada|huey', df_diatoms_in$`TRANSECT NAME`))
rows_delta_harnish_crescent_vg <- which(grepl('(?i)delta |harnish|crescent|von', df_diatoms_in$`TRANSECT NAME`))
rows_green_bowles <- which(grepl('(?i)green|bowles', df_diatoms_in$`TRANSECT NAME`))

# extract community data into wide format
df_comm_wide <- df_diatoms_in %>% 
  dplyr::select(-c(`record_no`,`region`,
                   `SAMPLE ID`, `ACCESSION NUMBER`, `TRANSECT NAME`,
                   `STREAM ABBREVIATION`, LATITUDE_DD, LONGITUDE_DD,                         
                   `COLLECTION DATE`, SEASON, `MAT TYPE`,                             
                   `AFDM (mg/cm2)`, `CHL-A (ug/cm2)`,
                   `TOTAL_ANNUAL_DISCHARGE (L)`, `VALVES COUNTED PER SAMPLE`, 
                   Unknown, `Marine diatom fragments`))

# fix species/column names so they don't break MCSim -- which requires r-friendsly col names
df_comm_wide <- df_comm_wide[,colSums(df_comm_wide) > 0]
names(df_comm_wide) <- names(df_comm_wide) %>%
  gsub(' ','_', .) %>%
  gsub('\\#','', .) %>%
  gsub('\\.','_', .) %>% tolower()
spp_list <- names(df_comm_wide)

# hellinger gransform data for OMI analysis
df_comm_wide_hell <- decostand(df_comm_wide[,spp_list], 'hellinger')

# extract sample info for observed records
df_sample_info <- df_diatoms_in %>% 
  dplyr::select(`SAMPLE ID`, `ACCESSION NUMBER`, `TRANSECT NAME`,
                `STREAM ABBREVIATION`, LATITUDE_DD, LONGITUDE_DD,                         
                `COLLECTION DATE`, SEASON, `MAT TYPE`,                             
                `AFDM (mg/cm2)`, `CHL-A (ug/cm2)`,
                `TOTAL_ANNUAL_DISCHARGE (L)`, `VALVES COUNTED PER SAMPLE`, 
                Unknown, `Marine diatom fragments`) 

# transform data and rename col names
df_sample_info <- df_sample_info %>% mutate(
  `MAT TYPE` = tolower(`MAT TYPE`),
  chl_a_log1p = log1p(`CHL-A (ug/cm2)`),
  afdm_log1p = log1p(`AFDM (mg/cm2)`),
  discharge_log1p = log1p(`TOTAL_ANNUAL_DISCHARGE (L)`),
  mat_color_orange = as.numeric(`MAT TYPE` == 'orange'),
  mat_color_black = as.numeric(`MAT TYPE` == 'black'))


# extract environmental variables to go into PCA
df_env <- df_sample_info %>% dplyr::select(
  chl_a_log1p, 
  # afdm_log1p, 
  # discharge_log1p,
  # mat_color_orange,
  # mat_color_green,
  # mat_color_red,
  mat_color_black)

# environmental PCA used to environmental axis/gradient
dudi_pca_env <- dudi.pca(df_env, 
                         scale = TRUE,
                         scan = FALSE,
                         nf = 1)

# calculate species habitat preferences and niche widths
niche_species <- niche(
  dudi_pca_env,
  Y = as.data.frame(df_comm_wide_hell),
  scann = FALSE
)

# extract niche position table from dudi object
df_niche <- data.frame(
  niche.pos = niche_species$li,
  as.data.frame(niche.param(niche_species))
)

#########################################
#########################################
#########################################
#########################################

# read in NDVI data
df_ndvi_in <- read_csv('Supplementary R Code/df_fryxell_basin_corrected_NDVI.csv') %>% slice(1:500) %>%
  arrange(desc(NDVI))

# rescale NDVI metrics to have distribution that is similar to chl a data
df_ndvi <- df_ndvi_in %>%
  mutate(NDVI_01 = (NDVI -  min(NDVI)) / (max(NDVI) - min(NDVI)),
         NDVI_rescaled = NDVI_01 * 150,
         NDVI_rescaledf_log1p = log1p(NDVI_rescaled))

# plot NDVI metrix on arbitrary xy map
df_ndvi %>% 
  ggplot(aes(row, col, 
             color = NDVI_rescaledf_log1p,
             size = NDVI_rescaledf_log1p)) +
  # geom_tile() +
  geom_point() +
  scale_color_gradient2(low = 'white', mid = scales::muted('blue'), high = 'green',
                        midpoint = 2.5)

# compare distributions for NDVI and chl a
par(mfrow = c(1,2))
hist(log1p(df_ndvi$NDVI_rescaled))
hist(df_env$chl_a_log1p)

#############################
# -- sample empirical data recrods to populate initial time step of simulation
df_simulation_site_info <- df_ndvi %>% mutate(
  obs_sample_row = dplyr::case_when(
    region == 'canada_huey' ~ sample(rows_canada_huey, nrow(df_ndvi), replace = TRUE),
    region == 'delta_harnish_crescent_vg' ~ sample(rows_delta_harnish_crescent_vg, nrow(df_ndvi), replace = TRUE),
    region == 'green_bowles' ~ sample(rows_green_bowles, nrow(df_ndvi), replace = TRUE)
  ),
  chl_a = df_sample_info$`CHL-A (ug/cm2)`[obs_sample_row],
  chl_a_log1p = df_sample_info$chl_a_log1p[obs_sample_row],
  transect_name = df_sample_info$`TRANSECT NAME`[obs_sample_row],
  mat_type = df_sample_info$`MAT TYPE`[obs_sample_row],
  afdm = df_sample_info$`AFDM (mg/cm2)`[obs_sample_row],
  Ef_scores = dudi_pca_env$li$Axis1[obs_sample_row])

#######################################
#######################################
#######################################
#######################################
# set up parallel

# n_cores <- 2
my_cl <- parallel::makeCluster(n_cores)

parallel::clusterEvalQ(
  my_cl,
  {
    options(stringsAsFactors = FALSE)
    library(MCSim)
    library(ade4)
    library(vegan)
    library(tidyverse)
    library(parallel)
  }
)

#######################################
# sim params
###################

design_matrix <- data.frame(
  niche_breadth_scaling_coef = sample(c(1, 4, 20), n_sims, replace = TRUE),
  m_immigration = sample(c(.001, .005, .009, .01, .05, .09, .1, .5, .9), n_sims, replace = TRUE), 
  m_immigration_black_mats = sample(c(.001, .005, .009, .01, .05, .09, .1, .5, .9), n_sims, replace = TRUE),
  nu = 10^(-1*sample(1:6, n_sims, replace = TRUE)),
  dispersal_kernel_slope = 10^sample(1:6, n_sims, replace = TRUE),
  i_sim = 1:n_sims) %>% distinct()

write_csv(design_matrix, 'design_matrix.csv')

################################
################################
################################
################################

# wrapper function to run simulations in parallel
# the default parameter settings here were for testing only and 
# are not used in the actual simulations (see the fxn call using 
# clusterMap below). 
fnx <- function(
  niche_breadth_scaling_coef = 4,
  m_immigration = .01, 
  m_immigration_black_mats = .05,  
  nu = .01,
  dispersal_kernel_slope = 200,
  i_sim = 1,
  ###### fixed ###########
  # df_ndvi = df_ndvi, #long form coordinages and ndvi values to use in landscape
  d_simulation_site_info = df_simulation_site_info,
  d_niche = df_niche,
  d_comm_wide = df_comm_wide,
  d_sample_info = df_sample_info,

  # group records by "region"
  rows_canada_huey = rows_canada_huey,
  rows_delta_harnish_crescent_vg = rows_delta_harnish_crescent_vg,
  rows_green_bowles = rows_green_bowles, 
  
  # other vars
  JL_max = 10000,
  JL_min = 50,
  invasion_limit = 10,
  n_timesteps = 100,
  path_to_read_data = getwd()){
  
  ##############################
  require(ade4)
  require(tidyverse)
  require(MCSim)
  require(vegan)
  
  #############################
  # -- sample empirical data recrods to populate initial time step of simulation
  
  m_values_tmp <- ifelse(
    d_simulation_site_info$mat_type == 'black',
    m_immigration_black_mats, m_immigration
  )
  
  d_simulation_site_info$m_values <- m_values_tmp
  
  #############################
  # create d_comm_t0 to initialize simulation
  d_comm_t0 <- d_comm_wide[d_simulation_site_info$obs_sample_row, ]
  d_comm_t0 <- d_comm_t0[,colSums(d_comm_t0) > 0]
  
  # extract habitat preferences and niche widths for species pool for simulation
  d_niche_t0 <- d_niche[names(d_comm_t0), ]
  
  ########################################
  # -- sim IDs
  scenario_description <- paste0(
    'METACOMM_nicheScaling-',niche_breadth_scaling_coef,
    '_m-',m_immigration,
    '_nu-', nu,
    '_w-',dispersal_kernel_slope,
    '_rep-',i_sim)

  try({
    i_sim_label <- paste0(scenario_description,'_', Sys.time()%>%format('%Y-%m-%d_%H%M%S'))#sim.ID
    
    simulation_landscape <- MCSim::fn.make.landscape(
      site.coords = d_simulation_site_info[,c('row','col')],
      Ef = d_simulation_site_info$Ef_scores,
      JL = round((d_simulation_site_info$NDVI_rescaled/max(d_simulation_site_info$NDVI_rescaled)) * JL_max) + JL_min,
      m = d_simulation_site_info$m_values,
      site.info = d_simulation_site_info)
    
    simoutput <- MCSim::fn.metaSIM(
      landscape = simulation_landscape,
      output.dir.path = paste0(path_to_read_data,'/SIM_RESULTS'),
      scenario.ID = 'streamMetacommunities',
      sim.ID = i_sim_label,
      trait.Ef = d_niche$Axis1,
      trait.Ef.sd = sqrt(d_niche$Tol) * niche_breadth_scaling_coef,
      J.t0 = as.data.frame(d_comm_t0),
      n.timestep = n_timesteps,
      W.r = dispersal_kernel_slope,
      nu = nu,
      speciation.limit = invasion_limit,
      save.sim = TRUE
    )
    
    print(paste('i_sim :',i_sim))
  })
} #END fnx

# ########################################
# ########################################
# ########################################
# ########################################


# test run commented out
# i <- 1
# fnx(niche_breadth_scaling_coef = design_matrix$niche_breadth_scaling_coef[i],
#     m_immigration = design_matrix$m_immigration[i],
#     m_immigration_black_mats = design_matrix$m_immigration_black_mats[i],
#     nu = design_matrix$nu[i],
#     dispersal_kernel_slope = design_matrix$dispersal_kernel_slope[i],
#     i_sim = design_matrix$i_sim[i],
#     path_to_read_data = getwd())


# use clusterMap to run simulations in parallel
clusterMap(
  cl = my_cl,
  fun = fnx,
  niche_breadth_scaling_coef = design_matrix$niche_breadth_scaling_coef,
  m_immigration = design_matrix$m_immigration,
  m_immigration_black_mats = design_matrix$m_immigration_black_mats,
  nu = design_matrix$nu,
  dispersal_kernel_slope = design_matrix$dispersal_kernel_slope,
  i_sim = design_matrix$i_sim,
  
  MoreArgs = list(
    ###### fixed ###########
    d_simulation_site_info = df_simulation_site_info,
    d_niche = df_niche,
    d_comm_wide = df_comm_wide,
    d_sample_info = df_sample_info,
    
    # group records by "region"
    rows_canada_huey = rows_canada_huey,
    rows_delta_harnish_crescent_vg = rows_delta_harnish_crescent_vg,
    rows_green_bowles = rows_green_bowles, 
    
    path_to_read_data = getwd()))


###########################
# stop parallel clusters
############
parallel::stopCluster(my_cl)
####################
