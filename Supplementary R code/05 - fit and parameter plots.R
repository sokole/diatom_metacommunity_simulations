rm(list = ls())
gc()

library(tidyverse)

#####################################################################
# get file names
file_name_list <- list.files()
file_name_list <- file_name_list[!grepl('q0', file_name_list)]

# find filename to load take most recent date if multiple
file_name_in <- sort(file_name_list[grepl('d_metacomm_aggregate_stats_', file_name_list)], decreasing = TRUE)[1]
d_agg <- readr::read_csv(file_name_in)


# find filename to load take most recent date if multiple
file_name_in <- sort(file_name_list[grepl('d_metacomm_fit_stats_aggregate_', file_name_list)], decreasing = TRUE)[1]
d_agg_fit <- read_csv(file_name_in)
# d_temporal <- read_csv('d_metacomm_stats_2018-06-10.csv')


# get metadata
d_sim_metadata <- read_csv('d_sim_metadata_streamMetacommunities.csv')
d_sim_metadata <- d_sim_metadata %>% 
  dplyr::rowwise() %>%
  mutate(
    m_max = max(m_black_mats, m_other),
    m_min = min(m_black_mats, m_other),
    m_ratio_abs = m_max/m_min) %>%
  mutate(
    `Metacomm. Dynamic` = case_when(
      niche_scaling == 1 ~ 'SSS',
      niche_scaling == 4 ~ 'MSS',
      niche_scaling == 20 ~ 'NM'
    ),
    fit_type = 'All')

d_sim_metadata$`Metacomm. Dynamic` <- d_sim_metadata$`Metacomm. Dynamic` %>%
  factor(levels = c('NM','MSS','SSS'), ordered = TRUE)

# join fit data
d_joined_metadata_aggfit <- d_sim_metadata %>%
  left_join(d_agg_fit)

# separate by dispersal type


  
d_strss <- d_joined_metadata_aggfit %>% filter(niche_scaling == 1)
d_modss <- d_joined_metadata_aggfit %>% filter(niche_scaling == 4)
d_neutral <- d_joined_metadata_aggfit %>% filter(niche_scaling == 20)

d_strss_best_fit <- d_strss %>%
  filter(chi_sq_tot <= quantile(d_strss$chi_sq_tot, probs = 0.05)) %>%
  mutate(fit_type = 'Top 5%')

d_modss_best_fit <- d_modss %>%
  filter(chi_sq_tot <= quantile(d_modss$chi_sq_tot, probs = 0.05)) %>%
  mutate(fit_type = 'Top 5%')

d_neutral_best_fit <- d_neutral %>%
  filter(chi_sq_tot <= quantile(d_neutral$chi_sq_tot, probs = 0.05)) %>%
  mutate(fit_type = 'Top 5%')

d_all <- d_joined_metadata_aggfit %>%
  bind_rows(d_strss_best_fit) %>%
  bind_rows(d_modss_best_fit) %>%
  bind_rows(d_neutral_best_fit)

d_plot <- d_all %>%
  rename(
    m_ratio = m_ratio_black_over_other,
    m_black = m_black_mats,
    m_orange = m_other) %>%
  gather('Parameter', 'Value', 
         chi_sq_tot,
         nu, w, 
         m_black, m_orange, m_ratio, m_ratio_abs) %>%
  select(scenario_id, sim_id, niche_scaling, `Metacomm. Dynamic`, fit_type,
         Parameter, Value)
  

###########################################################
###########################################################
# 
# # plot fit stat for all and best fit mods
# 
# FIG_boxplots_fit <- d_plot %>% 
#   filter(Parameter == 'chi_sq_tot') %>%
#   ggplot(
#     aes(
#       y = log(Value), 
#       x = fit_type,
#       color = fit_type)) +
#   geom_boxplot() +
#   facet_grid(. ~ `Metacomm. Dynamic`, scales = 'free') +
#   ylab('Deviation from initial condition\n(log Chi Sq. t_final - t_0)\n') +
#   xlab('Simulation scenarios') +
#   theme_bw() +
#   theme(legend.position="none")
# 
# # view plot
# graphics.off()
# windows(width = 4, height = 3, pointsize = 10)
# print(FIG_boxplots_fit)
# 
# # save plot 
# ggplot2::ggsave('FIG_boxplots_fit.jpg', FIG_boxplots_fit, 
#                 width = 4, height = 3, units = 'in')

#########
# plot fit as density

FIG_density_fit <- d_plot %>% 
  filter(Parameter == 'chi_sq_tot') %>%
  rename(`Sim. outcomes` = fit_type) %>%
  ggplot(
    aes(
      x = log(Value), 
      # x = fit_type,
      color = `Sim. outcomes`,
      fill = `Sim. outcomes`)) +
  geom_density(alpha = 0.25) +
  facet_grid(`Metacomm. Dynamic` ~ ., scales = 'free') +
  xlab(expression(paste('Deviation from initial condition ','(log ',Chi[sim]^2,')'))) +
  ylab('Density') +
  theme_bw()
# theme(legend.position="none")


# view plot
graphics.off()
windows(width = 4, height = 3.5, pointsize = 10)
print(FIG_density_fit)

# save plot 
ggplot2::ggsave('FIG_density_fit.jpg', FIG_density_fit, 
                width = 4, height = 3.5, units = 'in',
                dpi = 600)

###########################################################
###########################################################


# KS test -- is one distribution a random subset of the next?

d_distribution_comparison_stats <- data.frame()
for(i_meta_dyn in unique(d_plot$`Metacomm. Dynamic`)){
  for(i_param in unique(d_plot$Parameter)){
    
    d_i <- d_plot %>%
      dplyr::filter(
        `Metacomm. Dynamic` == i_meta_dyn, 
        Parameter == i_param)
    
    result_i <- ks.test(
      d_i %>% filter(fit_type == 'Top 5%') %>% select(Value) %>% unlist(use.names = FALSE),
      d_i %>% filter(fit_type == 'All') %>% select(Value) %>% unlist(use.names = FALSE))
    
    d_summary_i <- data.frame(
      `Metacomm. Dynamic` = i_meta_dyn, 
      Parameter = i_param,
      D = result_i$statistic,
      p_val = result_i$p.value,
      p_val_rounded = result_i$p.value %>% round(3)
    )
    
    d_distribution_comparison_stats <- d_distribution_comparison_stats %>%
      bind_rows(d_summary_i)
  }
}

write_csv(d_distribution_comparison_stats, 'RESULTS_d_distribution_comparison_stats_KS_test.csv')  



# plot parameter values for best and all models

FIG_boxplots_params <- d_plot %>% 
  rename(`Sim. outcomes` = fit_type) %>%
  filter(Parameter %in% c('nu','w','m_ratio')) %>%
  mutate(Parameter = gsub('m_ratio','m[black]/m[orange]',Parameter)) %>%
  ggplot(
    aes(
      x = log(Value), 
      # x = fit_type,
      color = `Sim. outcomes`,
      fill = `Sim. outcomes`)) +
  geom_density(alpha = 0.25,
               adjust = 1.5) +
  facet_grid(`Metacomm. Dynamic` ~ Parameter, scales = 'free',
             labeller = label_parsed) +
  ylab('Density') +
  xlab('Log parameter value') +
  theme_bw()
  # theme(legend.position="none")

# view plot
graphics.off()
windows(width = 6, height = 3.5, pointsize = 10)
print(FIG_boxplots_params)

# save plot 
ggplot2::ggsave('FIG_densplot_params.jpg', FIG_boxplots_params, 
                width = 6, height = 3.5, units = 'in',
                dpi = 600)
