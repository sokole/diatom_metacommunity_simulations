# diatom_metacommunity_simulations/Supplementary R Code

Supplementary R Code and intermediate files produced by these R scripts for:
Sokol ER, Barrett JE, Kohler TJ, McKnight DM, Salvatore MR and Stanish LF (2020) Evaluating Alternative Metacommunity Hypotheses for Diatoms in the McMurdo Dry Valleys Using Simulations and Remote Sensing Data. Front. Ecol. Evol. 8:521668. doi: 10.3389/fevo.2020.521668

R scripts:
 - 01 - create_ndvi_landscape.R - imports a processed remote image and calculates NDVI for areas of interest in the Fryxell Basin in Taylor Valley, Antarctica. 
 - 02 - parallel sims.R - imports NDVI data and MCM LTER archive data of diatoms counts in cyanobacteria mat cores sampled from stream mats in Fryxell Basin in Taylor Valley in Antarctica, conducts an OMI analysis, and then creates 10 000 metacommunity simulations using MCSim (v0.4.9). 
 - 03 - sim metadata.R - creates a metadata file with information about each simulation scenario parsed from its output file name. 
 - 04 - calc metacomm aggregate best fit stats.R - calculates alpha, beta, and gamma diversity for each simulation outcome at timesteps 1 and 100. Then calculates Chi-squared fit statistics as described in the article cited above.
 - 05 - fit and parameter plots.R - creates plots of (1) simulation Chi square statistics to identify the most plausible simulation scenarios and (2) a plot of parameter values across simulations.
 - 06 - summary plots.R - creates boxplots for alpha, beta, and gamma diversity faceted by simulation scenario at the end (timestep 100) of each simulation.
 
Intermediate data files produced by the above R scripts
 - d_metacomm_aggregate_stats_2019-10-04.csv
	- scenario_id - character, user defined name for a set of MCSim simulations
	- sim_id - character, unique simulation ID that includes parameter values and a time stamp
	- sim_filename - character, path to simulation output file
	- timestep - integer, timestep in the simulation for which the diversity metrics are being reported	
	- alpha_q1 - real, order q=1 alpha diveristy calculated for the metacommunity at the given timestep following Jost (2007)
	- beta_q1 - real, order q=1 beta diveristy calculated for the metacommunity at the given timestep following Jost (2007)	
	- gamma_q1 - real, order q=1 gamma diveristy calculated for the metacommunity at the given timestep following Jost (2007)
 - d_metacomm_fit_stats_aggregate_2019-10-04.csv
	- sim_id - character, unique simulation ID that includes parameter values and a time stamp
	- alpha_q1_length - integer, number of timesteps for which the given diversity metric was calculated for a given sim_id, should always be 2 because diversity metrics should only be calculated for timesteps 1 and 100. 
	- beta_q1_length - integer, number of timesteps for which the given diversity metric was calculated for a given sim_id, should always be 2 because diversity metrics should only be calculated for timesteps 1 and 100. 	
	- gamma_q1_length - integer, number of timesteps for which the given diversity metric was calculated for a given sim_id, should always be 2 because diversity metrics should only be calculated for timesteps 1 and 100. 	
	- alpha_q1_diff - real, change in given diveristy metric from timestep 1 to timestep 100	
	- beta_q1_diff - real, change in given diveristy metric from timestep 1 to timestep 100	
	- gamma_q1_diff - real, change in given diveristy metric from timestep 1 to timestep 100	
	- alpha_chi_sq - real, alpha compoment of the chi square metric, smaller values indicate less deviation from initial value during the course of the simulation	
	- beta_chi_sq - real, beta compoment of the chi square metric, smaller values indicate less deviation from initial value during the course of the simulation	
	- gamma_chi_sq - real, gamma compoment of the chi square metric, smaller values indicate less deviation from initial value during the course of the simulation	
	- chi_sq_tot - real, chi square metric to quantify deviation in alpha, beta, and gamma diversity at timestep 100 from initial condition for each simulation
 - d_sim_metadata_streamMetacommunities.csv
	- scenario_id - character, user defined name for a set of MCSim simulations
	- sim_file_name - character, path to simulation output file
	- sim_id - character, unique simulation ID that includes parameter values and a time stamp
	- niche_scaling	- real, niche scaling parameter
	- nu - real, invasion parameter	
	- w	- real, dispersal kernel slope
	- m_black_mats - real, immigration parameter for patches designated to represent black mats	
	- m_other - real, immigration parameter for patches that are not black mats (i.e., organe mats)	
	- proportion_black_mats	- real, proportion of patches (pixels) in the simulation designated to be black mats
	- proportion_orange_mats - real, proportion of patches (pixels) in the simulation designated to be orange mats	
	- proportion_red_mats - real, proportion of patches (pixels) in the simulation designated to be red mats (0 in this study)	
	- proportion_green_mats - real, proportion of patches (pixels) in the simulation designated to be green mats (0 in this study)	
	- proportion_biomass_black_mats - real, proportion of diatom biomass (i.e., individuals in the simulation) designated to reside in black mats	
	- proportion_biomass_orange_mats - real, proportion of diatom biomass (i.e., individuals in the simulation) designated to reside in orange mats		
	- proportion_biomass_red_mats - real, proportion of diatom biomass (i.e., individuals in the simulation) designated to reside in red mats (0 in this study)		
	- proportion_biomass_green_mats - real, proportion of diatom biomass (i.e., individuals in the simulation) designated to reside in green mats (0 in this study)		
	- m_ratio_black_over_other - real, ratio of m_black_mats / m_other
 - df_fryxell_basin_corrected_NDVI.csv - attributes of patches used in the metacommunity simulation, derived from a remote image of Fryxell Basin
	- row - integer, veritile position of pixel in raster image
	- col - integer, horizontal position of pixel in raster image
	- NDVI - real, Normalized Difference Vegetation Index (NDVI) calculated for pixels in the raster image, see main article text for details
	- NDVI_log1p - real, Transformed NDVI score, log(NDVI + 1)
	- NDVI_log1p_01 - real, rescaled NDVI_log1p to have min 0 and max 1
	- region - character, code identifies the region in Fryxell Basin in which the pixel is located, where the region identified as the name of the stream subcatchment(s) 