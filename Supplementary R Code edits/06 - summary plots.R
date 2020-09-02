#' R Script 06 for:
#' Sokol ER, Barrett JE, Kohler TJ, McKnight DM, Salvatore MR and Stanish LF (2020)
#'  Evaluating Alternative Metacommunity Hypotheses for Diatoms in the McMurdo Dry
#'  Valleys Using Simulations and Remote Sensing Data. Front. Ecol. Evol. 8:521668. 
#'  doi: 10.3389/fevo.2020.521668
#'
#' @author Eric R. Sokol \email{esokol@battelleecology.org}
#' 
#' This R script creates boxplots for alpha, beta, and gamma diversity faceted by 
#' simulation scenario at the end (timestep 100) of each simulation
 

# clean out workspace
rm(list = ls())
gc()


# load required packages
library(tidyverse)



# set paths for reading and writing data
# you may need to alter your path to the directory with this table
read_file_path <- "Supplementary R Code"
write_file_path <- "Supplementary R Code"


#####################################################################
# get file names
file_name_list <- list.files(read_file_path) %>% paste0(read_file_path,"/",.)


# find filename to load take most recent date if multiple
file_name_in <- sort(file_name_list[grepl("d_metacomm_aggregate_stats_", file_name_list)], decreasing = TRUE)[1]
d_agg <- readr::read_csv(file_name_in)


# find filename to load take most recent date if multiple
file_name_in <- sort(file_name_list[grepl("d_metacomm_fit_stats_aggregate_", file_name_list)], decreasing = TRUE)[1]
d_agg_fit <- read_csv(file_name_in)
# d_temporal <- read_csv("d_metacomm_stats_2018-06-10.csv")

d_sim_metadata <- read_csv(paste0(read_file_path,"/", "d_sim_metadata_streamMetacommunities.csv"))

d_sim_metadata <- d_sim_metadata %>% 
  dplyr::rowwise() %>%
  mutate(
    m_max = max(m_black_mats, m_other),
    m_min = min(m_black_mats, m_other),
    m_ratio_abs = m_max/m_min)



############################################################
d_agg_t1 <- d_agg %>% filter(timestep == min(timestep))

d_agg_t1_summary <- d_agg_t1 %>%
  summarize_at(c("alpha_q1", "beta_q1", "gamma_q1"), mean) %>%
  gather() %>% 
  mutate(
    Diversity = gsub("_q1","",key) %>% stringr::str_to_title()
  )

############################################################
d_plot <- left_join(d_sim_metadata,
                   d_agg %>% select(-scenario_id, -sim_filename),
                   by = "sim_id") %>%
  rename(Alpha = alpha_q1,
         Beta = beta_q1,
         Gamma = gamma_q1,
         m_ratio = m_ratio_black_over_other) %>%
  gather("Diversity", "Value", Alpha, Beta, Gamma) %>%
  mutate(Generation = as.factor(timestep),
         `Metacomm. Dynamic` = case_when(
           niche_scaling == 1 ~ "SSS",
           niche_scaling == 4 ~ "MSS",
           niche_scaling == 20 ~ "NM"
         ),
         `Dispersal Distance` = case_when(
           w == 1e6 ~ "Extremely Limited",
           w == 1e5 ~ "Limited",
           w == 100 ~ "Panmictic"
         )) %>% filter(Generation == 100) %>%
  mutate(
    Diversity = factor(Diversity, ordered = TRUE, 
                       levels = c("Gamma","Beta","Alpha")),
    `Metacomm. Dynamic` = factor(`Metacomm. Dynamic`, ordered = TRUE, 
                                 levels = c("SSS", "MSS", "NM"))
  )




# boxplots for alpha, beta, and gamma diversity faceted by simulation scenario
# at the end (timestep 100) of each simulation

FIG_agg_diversity_final_timestep <- d_plot %>% ggplot(
  # aes(log(m_ratio), Value, 
  aes(x = `Metacomm. Dynamic`,
      y = Value, 
      color = `Metacomm. Dynamic`)) +
  geom_boxplot(alpha = 0.1, size = .75) +
  # geom_smooth(method = "lm") +
  # geom_smooth() +
  # facet_grid(Diversity ~ `Metacomm. Dynamic`, scales = "free") +
  facet_grid(Diversity ~ ., scales = "free") +
  theme_bw() +
  ylab("Diversity (order q = 1 species equivalents)") +
  xlab("Metacomm. Dynamic") +
  # xlab(expression(paste("log(",m[black],"/",m[orange],")"))) +
  geom_hline(aes(yintercept = value), d_agg_t1_summary,
             linetype = "dashed",
             size = .75) +
  theme(legend.position="none")


# view plot
graphics.off()
windows(width = 3, height = 5, pointsize = 10)
print(FIG_agg_diversity_final_timestep)

# savePlot("FIG_agg_diversity_final_timestep", "pdf")
ggplot2::ggsave("FIG_agg_diversity_final_timestep.jpg", FIG_agg_diversity_final_timestep, 
                width = 4, height = 5,  units = "in",
                dpi = 600)
