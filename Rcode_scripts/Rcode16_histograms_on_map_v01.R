#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
# #remove everything in the working environment, without a warning!!
# rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(plyr) # Load plyr package
library(readr) # Load readr package
# # install package
# if(!require(rnaturalearth)){
#   install.packages("rnaturalearth")
#   install.packages("devtools")
#   #https://github.com/ropensci/rnaturalearthhires
#   remotes::install_github("ropensci/rnaturalearthhires")
#   install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")
# }

library("rnaturalearthhires")
# #define working directory
wd00 <- getwd()
#define input file  directory
wd01 <- "data/data_ddpcr_runs"
# paste together working directory and input file  directory
wd00_01 <- paste0(wd00,"/",wd01)

# make an output directory
wd16 <- "output16_ddPCR_and_qPCR_results_on_map_w_histograms"
# define the directory that holds the input file
wd06 <- "output06_compare_ddPCR_and_qPCR_results"

#make complete path to output dir
wd00_wd16 <- paste(wd00,"/",wd16,sep="")
#make complete path to input dir
wd00_wd06 <- paste(wd00,"/",wd06,sep="")
#Delete any previous versions of the output directory
unlink(wd00_wd16, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd16)

library(patchwork)

#
if(F){
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rworldmap)
  library(rworldxtra)
  library(sf)
  library(patchwork)
  library(ggrepel)  
  
  world <- rworldmap::getMap(resolution = "high")
  
  world <- sf::st_as_sf(world) %>%
    select(ADMIN)
  
}
# ----------------- options --------------------------

# set the map limits - these are also used to locate the inset plots
# they are assumed to be in lat/lon
map_bnds <- list(x0=-20, x1=42, y0=30, y1=75)
map_bnds <- list(x0=7, x1=17, y0=54.4, y1=58.2)
# set the projection / coordinate reference system to show the map in
map_crs <- 4326 # WGS84
#map_crs <- 3035

# should bar charts have an empty space for "missing" gene names 
show_zeros_in_bars <- TRUE

# name of output file 
png_name <- "plot.png"


# -----------------------------------------------------------
# evaluate on the projection 
if(map_crs==4326){
  plot_lat_lon <- TRUE
  plot_info <- "crs=4326 (WGS84) "
}else{
  plot_lat_lon <- FALSE 
  plot_info <- "crs=3035 (LAEA) "
}
map_crs <- sf::st_crs(map_crs) 

# if the map is plotted in a coordinate system other than WGS84, then we need
# to calculate the map corner points in this system

if(plot_lat_lon!=TRUE){
  map_crs <- sf::st_crs(map_crs)
  corners <- data.frame(x=c(map_bnds$x0, map_bnds$x0, map_bnds$x1, map_bnds$x1),
                        y=c(map_bnds$y0, map_bnds$y1, map_bnds$y1, map_bnds$y0))
  
  corners <- corners %>%
    sf::st_as_sf(coords=c("x","y"),
                 crs=sf::st_crs(4326)) 
  
  corners <- corners %>%
    sf::st_transform(crs=map_crs)
  
  corners <- sf::st_coordinates(corners) %>%
    as.data.frame()
  
  map_bnds$x0 <- min(corners$X)
  map_bnds$x1 <- max(corners$X)
  map_bnds$y0 <- min(corners$Y)
  map_bnds$y1 <- max(corners$Y)
}

# height and width of inset bar charts as fractions of the map height/width
inset_left <- 0.005 # horizontal offset from actual position of site (0-1)
inset_bottom <- 0.005 # vertical offset from actual position of site (0-1)
inset_height <- 0.12 # as fraction of height of plot (0-1)
inset_width <- 0.09 # as fraction of width of plot (0-1)

# ----------------- get sample data --------------------------
# # read sample data
df_c5 <- read.table(paste0(wd00_wd06,"/df_csvs05.2_v01.csv"),
                    header=T,sep=",",
                    fileEncoding = "latin1")
library(dplyr)
# rearrange the date frame to a longer format by gathering 
# the columns that starts with "l10"
# make a new column that holds the information 'dP' or 'qP',
# reflecting the 'l10.dP.cc' and the 'l10.qP.cc' column-header 
# that the value stems from.  This new have the header 'machine'.
df_c5 <- df_c5 %>%
  pivot_longer(cols = starts_with('l10.'),
               names_to = "machine",
               names_pattern = "l10\\.(.*)",
               values_to = "l10_copies")

# # exclude if it is NA
df_c5 <- df_c5[!is.na(df_c5$l10_copies),]
# get unique spc names
uniq_spcNm <- df_c5$spcNm %>% unique() %>% sort()
# get the number of species
Nspcs <- length(uniq_spcNm)
# In the ggplots , later on, a color is required for each species
# start by making a data frame
# get one range of colors
clr01 <- palette.colors(palette = "Okabe-Ito")
# get another range of colors
clr02 <- palette.colors(palette = "Polychrome")
# combine the two ranges of colors
clp <- c(clr01,clr02)
# add names to the palette colours so that they get matched to the spcNm values
names(clp) <- uniq_spcNm
# remove the "NA" value from the color palette
clp <- clp[!is.na(names(clp))]
# convert the gene names to factors
df_c5$spcNm <- factor(df_c5$spcNm, levels=uniq_spcNm)
# make sure the dplyr library is loaded
library(dplyr)
# count spcNm per site
spcNm_n <-  df_c5 %>%
  dplyr::group_by(smplNm, spcNm, dec_lat, dec_lon) %>%
  dplyr::summarise(n=n(), .groups="drop")
# get unique sites
sites <- spcNm_n %>%
  distinct(smplNm, dec_lat, dec_lon)
# remove previous versions of 'plot_info' if the object exists
if(exists(plot_info)==T){
  rm(plot_info)
}
# create a plot info string to use above the plot as a title
if(show_zeros_in_bars){
  # create a dataframe with 0-values where a gene name doesn't occur
  spcNm_n <- merge(sites,
                   distinct(spcNm_n,smplNm),
                   all=T) %>%
    left_join(spcNm_n, by=c("smplNm", "dec_lat", "dec_lon")) %>%
    mutate(n=ifelse(is.na(n),0,n))
  plot_info <- paste0(plot_info, "zeros shown in bars")
}else{
  plot_info <- paste0(plot_info, "zeros NOT shown in bars")
}

# if the 'l10_copies' is NA, then set it to 0
df_c5$l10_copies <- ifelse(is.na(df_c5$l10_copies), 0, df_c5$l10_copies)
# bar charts cannot be compared unless they all have the same y-axis
# to do this we need the maximum value possible across all sites
y_max <- max(spcNm_n$n) # 8
y_max <- max(df_c5$l10_copies) # 8
# ---- create main plot -------------------------------
sites <- sites[complete.cases(sites),]
# convert site positions to sf 
sites <- sites %>%
  sf::st_as_sf(coords=c("dec_lon", "dec_lat"), crs=sf::st_crs(4326))

if(plot_lat_lon!=TRUE){
  # we are not plotting in lat/lon
  # we need to get the site coordinates in the map coordinate system
  sites <- sites %>%
    sf::st_transform(crs=map_crs)
}
#sites <- sites[complete.cases(sites),]
coords <- sites %>%
  sf::st_coordinates() %>%
  as.data.frame()
names(coords) <- c("x","y")
sites <- sites %>%
  bind_cols(coords)
# load the library
library(ggplot2)
# plot the stations/sites
p_main <- ggplot() +
  geom_sf(data=world, colour=NA, fill="grey", alpha=0.5) + 
  geom_sf(data=sites, colour="black",shape=21,
          fill=alpha("white",0.1), size=2) + 
  coord_sf(xlim=c(map_bnds$x0,map_bnds$x1),
           ylim=c(map_bnds$y0,map_bnds$y1),
           crs=map_crs, datum = sf::st_crs(4326),
           expand = F)+
  theme_minimal() 

p <- p_main
# create a dummy data frame with points outside the map region
# this is used to create the legend but will not be visible
dummy <- data.frame(x=map_bnds$x1*1.1,y=map_bnds$y1*1.1)
dummy <- dummy %>%
  merge(data.frame(spcNm = factor(uniq_spcNm, levels=uniq_spcNm)),all=T)
# add the bar chart to the main plot
p <- p + 
  geom_bar(data=dummy, aes(x=x, fill=spcNm)) +
  scale_fill_manual(values=clp, name = "log10(copies/uL)\nfor assay") +
  theme(legend.position = "right") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + labs(subtitle=plot_info)
# see the plot
p
# loop through the study areas and add a bar chart to the main plot
for(site_name_i in unique(df_c5$smplNm)){
  print(site_name_i)
  #}
  # subset of species log10  copy numbers for this site
  dfi <- df_c5 %>%
    filter(smplNm==site_name_i)
  # plot the bar chart for the site
  pi <- ggplot() +
    geom_col(data=dfi, 
             aes(x=spcNm, 
                 y=l10_copies, fill=spcNm), na.rm=F) +
    scale_fill_manual(values=clp, guide="none") +
    scale_y_continuous(limits=c(0,y_max)) +
    coord_cartesian(expand=F) +
    theme_void() +
    theme(axis.line.x = element_line(linewidth=0.1, colour="lightgrey"))
  # find positions for the site as fraction of plot width/height
  sitei <- sites %>%
    filter(smplNm==site_name_i)
  
  x <- sitei$x[1]
  y <- sitei$y[1]
  
  # get the coordinates as fractions of distances from corner to corner
  x <- (x- map_bnds$x0) / (map_bnds$x1- map_bnds$x0)
  y <- (y- map_bnds$y0) / (map_bnds$y1- map_bnds$y0)
  
  
  # get the inset plot positions  
  x0 <- x + inset_left
  x1 <- x0 + inset_width * nrow(dfi)/length(uniq_spcNm)
  y0 <- y + inset_bottom
  y1 <- y0 + inset_height
  
  p <- p + 
    inset_element(pi, 
                  left = x0, 
                  bottom = y0,
                  right = x1,
                  top = y1)
  
}
# see the plot
p
# define the output file name
outfl <- paste0(wd00_wd16,"/Fig18_genes_used_in_geographical_regions_on_map.png")
#dev.off()
ggsave(p, filename=outfl, 
       height=210, 
       width = 297,
       units="mm", dpi=300)
#____________________________________________________________________________
#____________________________________________________________________________
# evaluate on the 'smpl.mnt' column, if it is lower or equal to 6, then it is
# 1st season, if it equal to 7 or higher, then it is 2nd season
df_c5$smpl.season <- ifelse(df_c5$smpl.mnt <= 6, "1st season", "2nd season")
unique(df_c5$smpl.season)
unique(df_c5$machine)
