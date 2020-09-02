#' R Script 01 for:
#' Sokol ER, Barrett JE, Kohler TJ, McKnight DM, Salvatore MR and Stanish LF (2020)
#'  Evaluating Alternative Metacommunity Hypotheses for Diatoms in the McMurdo Dry
#'  Valleys Using Simulations and Remote Sensing Data. Front. Ecol. Evol. 8:521668. 
#'  doi: 10.3389/fevo.2020.521668
#'
#' @author Eric R. Sokol \email{esokol@battelleecology.org}
#' 
#' This R script imports a processed remote image and calculates NDVI for areas
#' of interest in the Fryxell Basin in Taylor Valley, Antarctica. 


# required R packages
library(tidyverse)
library(raster)
library(ggplot2)


# Read in processed remote imagery data. File may be available upon request
# depending on data sharing permissions for past grantees who have used PGC 
# imagery with NSF OPP grants.  
image_tif <- raster('source_data/orthoWV02_15JAN211951582-M1BS-103001003ED2B400_u16ns3031_rad_atmcorr_refl_mask_NDVI_reduced.tif')
# plot(image_tif)

# str(image_tif)

# convert image to a matrix
image_txt <- as.matrix(image_tif)

# check number rows and cols in matrix produced from image
n_rows <- nrow(image_txt)
n_cols <- ncol(image_txt)

###############################################
# plot source image
######################
mat_to_plot <- image_txt

# replace negative values with 0
mat_to_plot[mat_to_plot<0] <- 0

# make color breaks
color_breaks <- c(0, quantile(mat_to_plot, c(0.75, 0.9, .99, .999), na.rm = TRUE))
color_list <- RColorBrewer::brewer.pal((length(color_breaks) - 1), 'RdYlGn')
color_list[1] <- 'light gray'

# view matrix as an image, save as a jpg
jpeg(filename="ndvi_full_scene.jpg")
image(mat_to_plot,
      col = color_list,
      breaks = color_breaks)
dev.off()
###############################################
# cropping out image areas that should not be used in the analysis
# e.g., glaciers, high elevation areas, areas outside of Fryxell Basin
image_tmp <- image_txt
x_thresh <- .525*n_rows
y_thresh <-  .525*n_cols
image_tmp[c(1:n_rows) > x_thresh, c(1:n_cols) > y_thresh] <- NA

x_thresh <- .75*n_rows
y_thresh <-  .35*n_cols
image_tmp[c(1:n_rows) > x_thresh, c(1:n_cols) > y_thresh] <- NA

y_max <- .75*n_cols
y_min <-  .35*n_cols
x_max <-  .75*n_rows
image_fryxell_basin <- image_tmp[c(1:n_rows) < x_max,
                                 c(1:n_cols) < y_max &  c(1:n_cols) > y_min]

###############################################
# plot cropped source image
######################
mat_to_plot <- image_fryxell_basin

# replace negative values with 0
mat_to_plot[mat_to_plot<0] <- 0

# make color breaks
color_breaks <- c(0, quantile(mat_to_plot, c(0.75, 0.9, .99, .999), na.rm = TRUE))
color_list <- RColorBrewer::brewer.pal((length(color_breaks) - 1), 'RdYlGn')
color_list[1] <- 'light gray'

# view cropped image and save as jpg
jpeg(filename="ndvi_cropped_fryxell_basin.jpg")
image(mat_to_plot,
      col = color_list,
      breaks = color_breaks)
dev.off()
###############################################

# only include top 1% NDVI values as pixels of interest
image_fryxell_basin_mask_low_vals <- image_fryxell_basin
image_fryxell_basin_mask_low_vals[
  image_fryxell_basin_mask_low_vals < quantile(image_fryxell_basin_mask_low_vals, .99, na.rm=TRUE)] <- NA

# plot image as jpg
# mat_to_plot <- image_txt
# mat_to_plot <- image_hoare_basin
# mat_to_plot <- image_fryxell_basin
# mat_to_plot <- image_outliers_masked
mat_to_plot <- image_fryxell_basin_mask_low_vals

color_breaks <- c(0, quantile(mat_to_plot, c(0.75, 0.9, .99, .999), na.rm = TRUE))
color_list <- RColorBrewer::brewer.pal((length(color_breaks) - 1), 'RdYlGn')
color_list[1] <- 'light gray'

# plot cropped with top 1% NDVI cells
jpeg(filename="ndvi_Fryxell_top_1_percent.jpg")
image(as.matrix(mat_to_plot),
      col = color_list,
      breaks = color_breaks)
dev.off()

###################################################################
# -- make raster matrix long, remove NAs, annotate records, and rescale NDVI metric
mat_in <- image_fryxell_basin_mask_low_vals
image_txt_long <- mat_in %>% as_tibble() %>%
  mutate(row = 1:nrow(mat_in)) %>%
  gather('col','NDVI', -row) %>% 
  mutate(col = as.integer(gsub('V','',col))) %>%
  filter(!is.nan(NDVI) & !is.na(NDVI)) %>%
  mutate(NDVI_log1p = log1p(NDVI),
         NDVI_log1p_01 = (NDVI_log1p - min(NDVI_log1p)) / ( max(NDVI_log1p) - min(NDVI_log1p)))

image_txt_long <- image_txt_long %>% mutate(
  region = case_when(
    row <= 1400 & col > 1250 ~ 'green_bowles',
    row <= 1750 & col %in% c(1001:1250) ~ 'delta_harnish_crescent_vg',
    row <= 2000 & col <= 1000 ~ 'delta_harnish_crescent_vg',
    row > 1400 & col > 1250 ~ 'canada_huey',
    row > 1700 & col %in% c(1001:1250) ~ 'canada_huey',
    row > 1700 & col %in% c(751:1000) ~ 'lost_seal_mcknight',
    row > 2200 & col %in% c(401:750) ~ 'lost_seal_mcknight',
    row > 2000 & col <= 400 ~ 'aiken',
    TRUE ~ 'delta_harnish_crescent_vg'))

# look at distributions of NDVI and rescaled NDVI values
hist(image_txt_long$NDVI)
hist(image_txt_long$NDVI_log1p_01)

# subsample records with top 500 NDVI values  
image_txt_long_top <- image_txt_long %>% arrange(desc(NDVI)) %>% slice(1:500)

# plot top 500 records
plot_top_500_NDVI_sites_by_region <- image_txt_long %>% 
  arrange(desc(NDVI)) %>% slice(1:500) %>%
  ggplot(aes(row, col, 
             color = region,
             size = NDVI)) +
  geom_point()

# save to pdf
pdf(file = 'ndvi_top_500_sites_by_region.pdf',
    width = 6, height = 3.5)
plot(plot_top_500_NDVI_sites_by_region)
dev.off()

# geom_tile() +
  # geom_point() +
  # scale_color_gradient2(low = 'white', mid = scales::muted('blue'), high = 'green',
  #                       midpoint = .225)

# write out csv of NDVI values
write_csv(image_txt_long_top, 'df_fryxell_basin_corrected_NDVI.csv')
