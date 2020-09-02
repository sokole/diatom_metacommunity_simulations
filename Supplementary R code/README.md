# diatom_metacommunity_simulations/Supplementary R Code

Supplementary R Code and intermediate files produced by these R scripts for:
Sokol ER, Barrett JE, Kohler TJ, McKnight DM, Salvatore MR and Stanish LF (2020) Evaluating Alternative Metacommunity Hypotheses for Diatoms in the McMurdo Dry Valleys Using Simulations and Remote Sensing Data. Front. Ecol. Evol. 8:521668. doi: 10.3389/fevo.2020.521668

R scripts:
 - 01 - create_ndvi_landscape.R
 - 02 - parallel sims.R
 - 03 - sim metadata.R
 - 04 - calc metacomm aggregate best fit stats.R
 - 05 - fit and parameter plots.R
 - 06 - summary plots.R
 
Intermediate data files produced by the above R scripts
 - d_metacomm_aggregate_stats_2019-10-04.csv
 - d_metacomm_fit_stats_aggregate_2019-10-04.csv
	- sim_id	
	- alpha_q1_length	
	- beta_q1_length	
	- gamma_q1_length	
	- alpha_q1_diff	
	- beta_q1_diff	
	- gamma_q1_diff	
	- alpha_chi_sq	
	- beta_chi_sq	
	- gamma_chi_sq	
	- chi_sq_tot

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
	- NDVI - real, Normalized Difference Vegetation Index calculated for pixels in the raster image, see main article text for details
	- NDVI_log1p - real, Transformed NDVI score, log(NDVI + 1)
	- NDVI_log1p_01 - real, rescaled NDVI_log1p to have min 0 and max 1
	- region - character, code identifies the region in Fryxell Basin in which the pixel is located, where the region identified as the name of the stream subcatchment(s) 
