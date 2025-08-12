#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(plyr)
library(readr)
library(rworldmap)
library(patchwork)
library(cowplot)
library(sf)

# Define working directory and input/output directories
wd00 <- getwd()
wd01 <- "data/data_ddpcr_runs"
wd00_01 <- paste0(wd00, "/", wd01)
wd16 <- "output16_ddPCR_and_qPCR_results_on_map_w_histograms"
wd06 <- "output06_compare_ddPCR_and_qPCR_results"
wd00_wd16 <- paste(wd00, "/", wd16, sep = "")
wd00_wd06 <- paste(wd00, "/", wd06, sep = "")
# unlink(wd00_wd16, recursive = TRUE)
# dir.create(wd00_wd16)

# Read and prepare data
df_c5 <- read.table(paste0(wd00_wd06, "/df_csvs05.2_v01.csv"),
                    header = T, sep = ",",
                    fileEncoding = "latin1")

df_c5 <- df_c5 %>%
  pivot_longer(cols = starts_with('l10.'),
               names_to = "machine",
               names_pattern = "l10\\.(.*)",
               values_to = "l10_copies")

df_c5 <- df_c5[!is.na(df_c5$l10_copies), ]
uniq_spcNm <- df_c5$spcNm %>% unique() %>% sort()
Nspcs <- length(uniq_spcNm)
clr01 <- palette.colors(palette = "Okabe-Ito")
#clr02 <- palette.colors(palette = "Polychrome")
clr02 <- palette.colors(palette = "Alphabet")
clp <- c(clr01, clr02)
names(clp) <- uniq_spcNm
clp <- clp[!is.na(names(clp))]
df_c5$spcNm <- factor(df_c5$spcNm, levels = uniq_spcNm)

# Evaluate on the 'smpl.mnt' column, if it is lower or equal to 6, then it is
# 1st season, if it is equal to 7 or higher, then it is 2nd season
df_c5$smpl.season <- ifelse(df_c5$smpl.mnt <= 6, "1st season", "2nd season")

spcNm_n <- df_c5 %>%
  dplyr::group_by(smplNm, spcNm, dec_lat, dec_lon) %>%
  dplyr::summarise(n = n(), .groups = "drop")

sites <- spcNm_n %>%
  distinct(smplNm, dec_lat, dec_lon)
# check if the plot_info variable exists, if it is remove it
if(exists("plot_info")==T){rm(plot_info)}

show_zeros_in_bars <- TRUE
if (show_zeros_in_bars) {
  spcNm_n <- merge(sites,
                   distinct(spcNm_n, smplNm),
                   all = T) %>%
    left_join(spcNm_n, by = c("smplNm", "dec_lat", "dec_lon")) %>%
    mutate(n = ifelse(is.na(n), 0, n))
}

df_c5$l10_copies <- ifelse(is.na(df_c5$l10_copies), 0, df_c5$l10_copies)
y_max <- max(df_c5$l10_copies)
ymax <- y_max
sites <- sites[complete.cases(sites), ]
sites <- sites %>%
  sf::st_as_sf(coords = c("dec_lon", "dec_lat"), crs = sf::st_crs(4326))

coords <- sites %>%
  sf::st_coordinates() %>%
  as.data.frame()
names(coords) <- c("x", "y")
sites <- sites %>%
  bind_cols(coords)
# if the world map is not loaded, load it
if(exists("world")==F){
  world <- rworldmap::getMap(resolution = "high")
  world <- sf::st_as_sf(world) %>%
    select(ADMIN)}
# Define map boundaries and coordinate reference system
map_bnds <- list(x0 = 7, x1 = 17, y0 = 54.4, y1 = 58.2)
map_crs <- sf::st_crs(4326)
# Get the global set of species names
global_uniq_spcNm <- unique(df_c5$spcNm)
# height and width of inset bar charts as fractions of the map height/width
inset_left <- 0.005 # horizontal offset from actual position of site (0-1)
inset_bottom <- 0.005 # vertical offset from actual position of site (0-1)
inset_height <- 0.52 # as fraction of height of plot (0-1)
inset_width <- 0.16 # as fraction of width of plot (0-1)


# should bar charts have an empty space for "missing" gene names 
show_zeros_in_bars <- TRUE
# data_subset <- df_c5.s1.m1
# mxuch <- mxuniqch
map_bnds_y0.rd <- floor(map_bnds$y0) 
map_bnds_y1.ru <- ceiling(map_bnds$y1)
# Function to create individual plots
create_plot <- function(data_subset,mxuch,y_max) {
  p_main <- ggplot() +
    geom_sf(data = world, colour = NA, fill = "grey", alpha = 0.5) +
    geom_sf(data = sites, colour = "black", shape = 21, fill = alpha("white", 0.1),
            size = 2) +
    coord_sf(xlim = c(map_bnds$x0, map_bnds$x1),
             ylim = c(map_bnds$y0, map_bnds$y1),
             crs = map_crs, datum = sf::st_crs(4326),
             expand = FALSE) +
    scale_y_continuous(breaks = seq(map_bnds_y0.rd, map_bnds_y1.ru, by = 1.0)) +
    #adjust the font size of the tick labels
    theme(
      axis.text.x = element_text(size = 5.2),  # Adjust the size of x-axis tick labels
      axis.text.y = element_text(size = 5.2)   # Adjust the size of y-axis tick labels
    ) +
    theme_minimal()
  
  p <- p_main
  # Calculate the maximum number of unique categories across all subplots
  max_unique_categories <- mxuch
  # Modify the loop to dynamically adjust the inset_width
  for (site_name_i in unique(data_subset$smplNm)) {
    print(site_name_i)
    #}
    #data_subset <- df_c5
    # Subset of gene counts for this site
    dfi <- data_subset %>%
      filter(smplNm == site_name_i)
    
    # Ensure the species factor uses the full set of species from the global data
    dfi$spcNm <- factor(dfi$spcNm, levels = global_uniq_spcNm)
    dfi$charNm <- dfi$spcNm 
    dfi$n <- dfi$l10_copies
    
    # Plot the bar chart for the site
    pi <- ggplot() +
      geom_col(data=dfi, aes(x = spcNm,
                             y = l10_copies,
                             fill = spcNm),
               width = 1,
               na.rm=F) +
      scale_fill_manual(values=alpha(clp,0.76), guide="none") +
      scale_y_continuous(limits=c(0,y_max)) +
      coord_cartesian(expand=F) +
      theme_void() +
      theme(axis.line.x = element_line(linewidth=0.1, colour="lightgrey"))
    
    # Find positions for the site as fraction of plot width/height
    sitei <- sites %>%
      filter(smplNm == site_name_i)
    
    x <- sitei$x[1]
    y <- sitei$y[1]
    
    # Get the coordinates as fractions of distances from corner to corner
    x <- (x - map_bnds$x0) / (map_bnds$x1 - map_bnds$x0)
    y <- (y - map_bnds$y0) / (map_bnds$y1 - map_bnds$y0)
    
    # Dynamically adjust the inset_width based on the number of unique categories in this subplot
    dynamic_inset_width <- inset_width * (length(unique(dfi$charNm)) / max_unique_categories)
    
    # Get the inset plot positions
    x0 <- x + inset_left
    x1 <- x0 + dynamic_inset_width*2
    y0 <- y + inset_bottom
    y1 <- y0 + inset_height
    
    p <- p +
      inset_element(pi,
                    left = x0,
                    bottom = y0,
                    right = x1,
                    top = y1)
  }
  
  
  return(p)
}
# Get unique machine values
unique_machines <- unique(df_c5$machine)
mxuniqch <- length(df_c5$spcNm %>% unique())
# subset the data frame by machine and season
df_c5.s1.m1 <- filter(df_c5, smpl.season == "1st season", machine == unique_machines[1])
df_c5.s1.m2 <- filter(df_c5, smpl.season == "1st season", machine == unique_machines[2])
df_c5.s2.m1 <- filter(df_c5, smpl.season == "2nd season", machine == unique_machines[1])
df_c5.s2.m2 <- filter(df_c5, smpl.season == "2nd season", machine == unique_machines[2])
#df_c5$l10_copies
# Create the four plots
plot1 <- create_plot(df_c5.s1.m1, mxuniqch,ymax)
plot2 <- create_plot(df_c5.s1.m2, mxuniqch,ymax)
plot3 <- create_plot(df_c5.s2.m1, mxuniqch,ymax)
plot4 <- create_plot(df_c5.s2.m2, mxuniqch,ymax)

# # Combine the four plots into a single plot
# combined_plot <- plot_grid(
#   #plot1 + theme(legend.position = "none"),
#   plot1 + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), units = "pt") ) ,
#   #theme(plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))
#   plot2 + theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), units = "pt") ),
#   plot3 + theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), units = "pt") ),
#   plot4 + theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), units = "pt") ),
#   #labels = "AUTO"
#   labels = c(letters[1:4]),
#   ncol = 2
# )

# Create a dummy plot to extract the legend
dummy_plot <- ggplot() +
  geom_col(aes(x = "", y = 0, fill = factor(uniq_spcNm))) +
  scale_fill_manual(values = clp, name = "log10(copies/uL)\nfor assay") +
  theme_void()

legend <- cowplot::get_legend(dummy_plot)

cmb_plt <- plot_grid(
  plot1 + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), units = "pt") ) ,
  plot2 + theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), units = "pt") ),
  plot3 + theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), units = "pt") ),
  plot4 + theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), units = "pt") ),
  ncol = 2,  # Organizing into a 2x2 grid
  labels = c("a", "b", "c", "d"),  # Adding subplot labels
  label_size = 14,
  rel_widths = c(1, 1),  # Reduce space for subplots (closer together)
  rel_heights = c(1, 1), # Reduce height for subplots (closer together)
  align = 'v',  # Align vertically to minimize the gap
  axis = "tb",  # Align axes at top and bottom
  hjust = -1,   # Move the labels to a better position
  vjust = 0.98,  # Adjust label position to move it up
  nrow = 2      # Force to two rows
  
)

combined_plot <- plot_grid(cmb_plt, 
                           legend, 
                           ncol = 2, 
                           rel_widths = c(1.7, 0.3),
                           rel_heights = c(2, 0.3),
                           align = "v")
# Save the combined plot
outfl <- paste0(wd00_wd16, "/Fig19_combined_plots_08.png")
ggsave(combined_plot, filename = outfl,
       height = 297 * 0.42,
       width = 210,
       units = "mm", dpi = 300)
