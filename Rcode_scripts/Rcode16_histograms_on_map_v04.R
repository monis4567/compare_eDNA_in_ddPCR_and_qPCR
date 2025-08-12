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
library(grid)

# Define working directory and input/output directories
wd00 <- getwd()
wd01 <- "data/data_ddpcr_runs"
wd00_01 <- paste0(wd00, "/", wd01)
wd16 <- "output16_ddPCR_and_qPCR_results_on_map_w_histograms"
wd06 <- "output06_compare_ddPCR_and_qPCR_results"
wd00_wd16 <- paste(wd00, "/", wd16, sep = "")
wd00_wd06 <- paste(wd00, "/", wd06, sep = "")
#Create a directory to obtain the data frames from
wd00_wd02 <- paste0(wd00,
                    paste0("/output13_",#,no_sbs,
                           "_my_df_modified_to_match",
                           "_G_Guri_setup"
                    ))
# unlink(wd00_wd16, recursive = TRUE)
# dir.create(wd00_wd16)
#define data directory
wd00.1 =paste0(wd00,"/data/MONIS6_2021_data")
# read in an xlsx file with all primer assays listed 
# for each species and the abbreviations
inf04 <- "list_of_specific_assays_MONIS6.xlsx"
wd00.1_inf04 <- paste0(wd00.1,"/",inf04 )
#library(xlsx)
df_dtc_asss <- openxlsx::read.xlsx(wd00.1_inf04,1)
# use the data frame with the detection assays to get the full
# length Latin Genus and Species names
df_dtc_asss$GnNm.l1 <- substr(df_dtc_asss$Genus, 1, 1)
df_dtc_asss$shrtNm <- paste0(df_dtc_asss$GnNm.l1,". ",df_dtc_asss$Species)
# Use plyr to select the columns 'shrtNm', 'Phyla', 'Class', 'AbbrvNm', 'Species'
# 'Genus', 'Organismegruppetype', 'Latinsk_navn', 'Navn_dansk' from the 'df_dtc_asss'
# data frame
df_dta <- dplyr::select(df_dtc_asss, shrtNm, Phyla, Class, AbbrvNm, Species,
                        Genus, Organismegruppetype, Latinsk_navn, Navn_dansk)
# limit the 'df_dta' data frame to remove rows where there is duplicate
# occurences of 'AbbrvNm'
df_dta <- df_dta[!duplicated(df_dta$AbbrvNm),]
#read in the tsv tables
df_ddpcr03 <- read.table(paste0(wd00_wd02,"/df_ddpcr03.tsv"),
                         sep = "\t",
                         header = T) 
df_qpcr03 <- read.table(paste0(wd00_wd02,"/df_qpcr03.tsv"),
                        sep = "\t",
                        header = T)
#in the 'df_qpcr03$Species' and the 'df_ddpcr03$Species' replace "Cragig"
# with 'Maggig', 
df_qpcr03$Species <- gsub("Cragig", "Maggig", df_qpcr03$Species)
df_ddpcr03$Species <- gsub("Cragig", "Maggig", df_ddpcr03$Species)
# get the unique species names from the 'df_ddpcr03' data frame
uqp_spc<- unique(df_qpcr03$Species)
udp_spc<- unique(df_ddpcr03$Species)
# combine the two unique species names in one vector
spcabbr.unq <- unique(c(uqp_spc,udp_spc))
# order them alphabetically
spcabbr.unq <- sort(spcabbr.unq)
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

#in the 'df_c5$spcNm' replace "Cragig"
# with 'Maggig', 
df_c5$spcNm <- gsub("Cragig", "Maggig", df_c5$spcNm)
# then use 'left_join' to merge the 'df_dta' data frame with the 'df_c5' using
# the 'AbbrvNm' in the 'df_dta' data frame and the 'spcNm' in the 'df_c5' data frame
# as the common column
df_c5 <- dplyr::left_join(df_c5, df_dta, by = c("spcNm" = "AbbrvNm"))


uniq_spcNm <- df_c5$spcNm %>% unique() %>% sort()
Nspcs <- length(uniq_spcNm)
clr01 <- palette.colors(palette = "Okabe-Ito")
clr02 <- palette.colors(palette = "Alphabet")
clp <- c(clr01, clr02)
names(clp) <- uniq_spcNm
clp <- clp[!is.na(names(clp))]
df_c5$spcNm <- factor(df_c5$spcNm, levels = uniq_spcNm)
# Create a color palette for species
clr03 <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#999999", "#AA0DFE", "#3283FE", "#85660D",
  "#782AB6", "#565656", "#1C8356", "#16FF32",
  "#F7E1A0"
)
# use the combined range of colors, to match up with the
# number of species in the data 'tibl_podd' data frame
clr04 <- clr03[1:length(unique(spcabbr.unq))]
col.f_spc <- clr04
spcNm <- spcabbr.unq
# combine in to data frame
df_clr <- as.data.frame(cbind(spcNm,col.f_spc))
# select only the 'shrtNm' and the 'Species' using filter
# and  keep only unique rows using 'distinct'
df_c5.sNm <- df_c5 %>% dplyr::select(shrtNm,spcNm) %>% distinct()
# use 'left_join' to merge the 'df_c5.sNm' data frame with the
# 'df_clr' data frame
df_clr <- dplyr::left_join(df_clr,df_c5.sNm, by = c( "spcNm"="spcNm" ))
#View(df_c5)
uniq_shrtNm <- df_c5$shrtNm %>% unique() %>% sort()


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
inset_height <- 0.44 # as fraction of height of plot (0-1)
inset_width <- 0.16 # as fraction of width of plot (0-1)

# should bar charts have an empty space for "missing" gene names
show_zeros_in_bars <- TRUE
# data_subset <- df_c5.s1.m1
# mxuch <- mxuniqch
map_bnds_y0.rd <- floor(map_bnds$y0)
map_bnds_y1.ru <- ceiling(map_bnds$y1)
# Function to create individual plots
create_plot <- function(data_subset, mxuch, y_max) {
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
      scale_fill_manual(values=alpha(setNames(df_clr$col.f_spc, df_clr$spcNm),0.76), guide="none") +
      
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
str(df_c5)
# subset the data frame by machine and season
df_c5.s1.m1 <- filter(df_c5, smpl.season == "1st season", machine == unique_machines[1])
df_c5.s1.m2 <- filter(df_c5, smpl.season == "1st season", machine == unique_machines[2])
df_c5.s2.m1 <- filter(df_c5, smpl.season == "2nd season", machine == unique_machines[1])
df_c5.s2.m2 <- filter(df_c5, smpl.season == "2nd season", machine == unique_machines[2])
# Create the four plots
plot1 <- create_plot(df_c5.s1.m1, mxuniqch, ymax)
plot2 <- create_plot(df_c5.s1.m2, mxuniqch, ymax)
plot3 <- create_plot(df_c5.s2.m1, mxuniqch, ymax)
plot4 <- create_plot(df_c5.s2.m2, mxuniqch, ymax)
#colnames(df_c5)
uniq_shrtNm <- df_c5$shrtNm %>% unique() %>% sort()
# define output directory and write out the table w colors for species
wdd <- "data/MONIS6_2021_data"
wd00_wdd <- paste0(wd00,"/",wdd)
fl.out <- "table_w_colors_for_spc.csv"
fl.out <- paste0(wd00_wdd,"/",fl.out)
write.table(df_clr, file=fl.out, sep=";",
            row.names = T)
# Create a dummy plot to extract the legend for species
dummy_plot_species <- ggplot() +
  geom_col(aes(x = "", y = 0, 
               fill = factor(uniq_shrtNm))) +
  scale_fill_manual(values = setNames(df_clr$col.f_spc, 
                    df_clr$shrtNm), name = "Species") +
  theme_void() +
  theme(legend.key.width = unit(0.4, 'cm'),  # Adjust the key width
        legend.key.height = unit(0.4, 'cm'), # Adjust the key height
        legend.text = element_text(size = 10,face = "italic")) # Adjust the text size
# Extract the legend for species
legend_species <- cowplot::get_legend(dummy_plot_species)
# Create a data frame for the height legend
height_df <- data.frame(height = seq(0, y_max, length.out = 100))
# Create a plot to represent the bar heights
height_plot <- ggplot(height_df, aes(x = 1, y = height, fill = height)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "blue"),
                       name = "log10(copies/uL)") +
  scale_y_continuous(breaks = seq(0, y_max, by = 2)) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.4, 'cm'))

# Extract the legend for height
legend_height <- cowplot::get_legend(height_plot)

# Combine the four plots into a single plot
cmb_plt <- plot_grid(
  plot1 + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), units = "pt")),
  plot2 + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), units = "pt")),
  plot3 + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), units = "pt")),
  plot4 + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), units = "pt")),
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

# Combine the plots and legends
combined_plot <- plot_grid(
  cmb_plt,
  plot_grid(legend_species, legend_height, ncol = 1, rel_heights = c(0.6, 0.4)),
  ncol = 2,
  rel_widths = c(1.7, 0.3),
  rel_heights = c(2, 0.3),
  align = "v"
)

# Save the combined plot
outfl <- paste0(wd00_wd16, "/Fig19_combined_plots_15.png")
ggsave(outfl, plot = combined_plot, height = 297 * 0.42, width = 210, units = "mm", dpi = 300)
