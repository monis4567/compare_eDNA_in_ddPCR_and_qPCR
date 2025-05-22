#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

library(dplyr)
library(tidyr)

#____________________________________________________________________________#
# R-code provided for the project:
# “comparison of qPCR and ddPCR”
# Authors: Steen Wilhelm Knudsen.
infile <- "table_w_ddpcr_results_v02.csv"
wd00   <- getwd()
#define input file  directory
wd01 <- "data/data_ddpcr_runs"
# paste together working directory and input file  directory
wd00_01 <- paste0(wd00,"/",wd01)
# make an output directory
wd12 <- "output12_map_iNat_records_and_ddpcr"
# make an output directory
wd13 <- "output13_heatmap_records_ddpcr_and_qpcr"
#make complete path to output dir
wd00_wd12 <- paste(wd00,"/",wd12,sep="")
pth_infl <- paste0(wd00_wd12,"/",infile)
# read in the data file
#df_A08 <- read.csv(pth_infl)
#make complete path to output dir
wd00_wd13 <- paste(wd00,"/",wd13,sep="")
#Delete any previous versions of the output directory
unlink(wd00_wd13, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd13)

#______________________________

#################################################################################

# define folder with previously stored csv file
wdout11 <- "output11_compare_all_ddpcr_and_qpcr"
# read in the data frame as a csv file
df_e12 <- read.csv(file=paste0(wd00,"/",wdout11,
                     "/table_w_ddpcr_and_qpcr_data_and_lat_lon_v03.csv"))
# make an empty column that can have information on how many
# uL was used as template in the ddPCR and the qPCR setup
df_e12$voltmpl_uL <- NA
# all ddPCR setups used 5 uL
df_e12$voltmpl_uL[(df_e12$mch=="ddPCR")] <- 5
# all qPCR setups used 3 uL
df_e12$voltmpl_uL[(df_e12$mch=="qPCR")] <- 3

#df_e12$Elueringvolumen.AE.buffer.uL # is the volume of elution buffer used for the extraction
#df_e12$voltmpl_uL # is the volume of template used
#df_e12$molcnt # is copies per uL counted
#df_e12$Vwf_mL # is the volume of water filtered in mL
# if the 'Elueringvolumen.AE.buffer.uL' represents the entire
# extraction, then the 'voltmpl_uL' used in the PCR is the fraction 
# of this entire extraction
frcextr <- df_e12$Elueringvolumen.AE.buffer.uL/df_e12$voltmpl_uL
# this 'frcextr' can then be multiplied with the 'molcnt' to get the
# total number of copies present in the extraction
cpinextr <-frcextr*df_e12$molcnt
# the 'cpinextr' can then be divided by the volume of water filtered
# to get the concentration of copies per mL
df_e12$cp_mL <- cpinextr/df_e12$Vwf_mL
# the 'cp_mL' can then be multiplied by 1000 to get the concentration of
# copies per L
df_e12$cp_L <- df_e12$cp_mL*1000

####################################################################################
# Start Appendix E
####################################################################################


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Use ipdw package to interpolate between marine sampling locations
# interpolate between sampling locations using coastlines as barriers
#-  as the fish swims, not as the crow flies!
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#code is prepared in 2024-Dec by Steen W. Knudsen
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#remove everything in the working environment, without a warning!!
#rm(list=ls())

#libr.path <- "/home/sknu003/uoa00029_runs/Rplot_tryout"
#.libPaths( c( libr.path, .libPaths()) )

#libr.path <- "/home/sknu003/uoa00029_runs/Rplot_tryout"
#libr.path <- "/scale_wlg_persistent/filesets/home/sknu003/R/x86_64-pc-linux-gnu-library/3.5"

#libr.path <- "/scale_wlg_persistent/filesets/home/sknu003/R/x86_64-pc-linux-gnu-library/3.6"
#.libPaths( c( .libPaths(), libr.path) )

#.libPaths()

#.libPaths( c( libr.path , .libPaths() ) )
#.libPaths()
#.libPaths(libr.path)
#.libPaths()
#chooseCRANmirror(graphics=FALSE)
#chooseCRANmirror(4)
#'chooseCRANmirror(graphics=FALSE, ind=4)'
#________________________________________________________________________________
#01 - use the two spatial dataframes as in this example    https://jsta.github.io/ipdw/articles/ipdw2.html

# also check out this website: https://globalfishingwatch.org/data-blog/working-with-our-downloadable-public-data-in-r/
# and this website: https://www.molecularecologist.com/2015/07/marmap/
#________________________________________________________________________________
# get the rgeos package
# if(!require(rgeos)){
#   install.packages("rgeos", repos='http://cran.us.r-project.org')
#   
# }
# remotes::install_version("rgeos", version = "0.6-4")
# library(rgeos)

# get the ipdw package
# if(!require(ipdw)){
#   install.packages("ipdw", repos='http://cran.us.r-project.org')
# }
library(ipdw)
# get the scales package
# if (!requireNamespace("scales", quietly=TRUE))
#   install.packages("scales", repos='http://cran.us.r-project.org')
library(scales)
# get the sf package
# if (!requireNamespace("sf", quietly=TRUE))
#   install.packages("sf", repos='http://cran.us.r-project.org')
library(sf)

# get the rnaturalearth package
# if (!requireNamespace("rnaturalearth", quietly=TRUE))
#   install.packages("rnaturalearth", repos='http://cran.us.r-project.org')
library(rnaturalearth)

#Read in the rgdal library
# if (!requireNamespace("rgdal", quietly=TRUE))
#   install.packages("rgdal", repos='http://cran.us.r-project.org')
#library(rgdal)

# get the sp package
# if (!requireNamespace("sp", quietly=TRUE))
#   install.packages("sp", repos='http://cran.us.r-project.org')
library(sp)

#https://www.rdocumentation.org/packages/biogeo/versions/1.0/topics/dms2dd
# get biogeo package to be able to use 'dms2dd' function
# if(!require(biogeo)){
#   install.packages("biogeo", repos='http://cran.us.r-project.org')
# }
# library(biogeo)
#https://www.rdocumentation.org/packages/rgdal/versions/1.3-6/topics/readOGR

# if(!require(spatstat)){
#   install.packages("spatstat", repos='http://cran.us.r-project.org')
# }
library(spatstat)
library(ggplot2)

# Read the .csv file with collection localities
# the lon-lat 02 positions are positions close to the sampling site, 
# but might not be
# positions in the 'sea' when the 'worldmap' in R is used. 
# This might be a problem
# for positions in fiords and narrow straits. 
# To be able to use these sampling locations
# an extra position further away in the sea is included. 
# This is the lon-lat 03 positions.
# The idea is that these positions further away
# will enable the R-package 'ipdw' to
# interpolate between the sampling locations, even though some of 
# the sampling locations
# might be in narrow straits and fiords.

# get the month abbreviation
df_e12$Dato_inds.mm2 <- month.abb[df_e12$Dato_inds.mm]

#count the number of season to loop over
no.of.seasons <- length(unique(df_e12$ssn.smpl))
# make a sequence of numbers to use in a data frame
no_for_season <- seq(1:no.of.seasons)
#get the names of the seasons -  to use in the loop below
categories.of.seasons <- sort(unique(df_e12$ssn.smpl))
# make names for the seasons
names.of.seasons <- c("spring","fall")
# bind to a data frame
seaons_nms_df <- as.data.frame(cbind(no_for_season,categories.of.seasons,names.of.seasons))
# make one of the columns numeric
seaons_nms_df$no_for_season <- as.numeric(seaons_nms_df$no_for_season)

#https://stackoverflow.com/questions/47418127/r-how-to-aggregate-some-columns-while-keeping-other-columns
#paste columns together
df_e12$spcAbbr.yr.smpln.ssn <- paste(df_e12$speciesabbr,
                               df_e12$Dato_inds.yy,
                               df_e12$smplNm,
                               df_e12$ssn.smpl,
                                sep=".")

###########################################################################
# Get the highest levels of eDNA
#get max value per group - this will be needed for
# the ipdw plots, where you want to
#set a max value on the legend for the heatmap
#https://stackoverflow.com/questions/25314336/extract-the-maximum-value-within-each-group-in-a-dataframe
df_mx_mlcnt <- aggregate(df_e12$molcnt, 
                                       by = 
                                         list(
                                           df_e12$Dato_inds.yy,
                                           df_e12$mch,
                                           df_e12$ssn.smpl,
                                           df_e12$speciesabbr), 
                                       max)
# also calculate the max value for the copy number per L
df_mx_cp_pL <- aggregate(df_e12$cp_L, 
                         by = 
                           list(
                             df_e12$Dato_inds.yy,
                             df_e12$mch,
                             df_e12$ssn.smpl,
                             df_e12$speciesabbr), 
                         max)

# alter the column headers
colnames(df_mx_mlcnt) <- c("Dato_inds.yy",
                           "mch",
                            "ssn.smpl",
                            "speciesabbr",
                           "mxmolcnt")

# alter the column headers
colnames(df_mx_cp_pL) <- c("Dato_inds.yy",
                           "mch",
                           "ssn.smpl",
                           "speciesabbr",
                           "mxcp_L")
# make all columns in the dataframe character
library(dplyr)
df_mx_mlcnt <- df_mx_mlcnt %>%
  mutate(across(everything(), as.character))
# do the same character change for all columns in the dataframe with
# maximum number of copies per L
df_mx_cp_pL <- df_mx_cp_pL %>%
  mutate(across(everything(), as.character))
# make the sampling year a character
df_e12$Dato_inds.yy <- as.character(df_e12$Dato_inds.yy)
#match back to original dataframe
df_e12 <- dplyr::left_join(df_e12,
                df_mx_mlcnt, 
                by = c("Dato_inds.yy","mch","ssn.smpl","speciesabbr"))
#match back to original dataframe , to also get the max number of copies per L
# per year, per machine per season and per species
df_e12 <- dplyr::left_join(df_e12,
                df_mx_cp_pL, 
                by = c("Dato_inds.yy","mch","ssn.smpl","speciesabbr"))
# make the maximum number of copies per L a numeric
df_e12$mxcp_L <- as.numeric(df_e12$mxcp_L)
df_e12$mxmolcnt <- as.numeric(df_e12$mxmolcnt)
# Get logarhitmic values of the maximum number of copies per L
df_e12$l10cp_L <- log10((df_e12$cp_L+1))
df_e12$l10mxcp_L <- log10((df_e12$mxcp_L+1))
###########################################################################
#try and get the coastline that is downloaded in a supporting
# directory
coastline10 <- ne_download(scale = 10, 
                       type = 'land', 
                       category = 'physical', 
                       destdir = paste(wd00_wd13,sep=""))
#check if the object with the coastline does not exists
if (!exists("coastline10"))
{
  #get the coastline
  # a scale close to '10' gives a fine detailed coast, but takes longer to calculate
  # a scale above to '60' gives a non-detailed coast, but is faster to calculate
  #coastline110 <- ne_download(scale = 110, type = 'land', category = 'physical')
  #coastline10 <- ne_download(scale = 110, type = 'land', category = 'physical')
  #coastline10 <- ne_load(scale = 10, type = 'land', category = 'physical', destdir = paste(wd00,wd09,sep=""))
  coastline10 <- ne_download(scale = 10, type = 'land', category = 'physical')
  #coastline10 <- ne_download(scale = 40, type = 'land', category = 'physical')
  
  #close the if test above - i.e. if the 'coastline10' object does not exist, then get it
}
###########################################################################

# Assigning CRS
#Note the CRS is different from the UTM CRS prepared above for this CRS based on lonlat
r2 <- sp::CRS("+init=epsg:4326") # the 4326 works for Northern Europe
# create crs object
epsg4326nCRS2 <- crs(r2)

# read in an xlsx file with all primer assays listed 
# for each species anmd the abbreviations
inf04 <- "list_of_specific_assays_MONIS6.xlsx"
# copy the path to the data folder
wd00_wddata <- paste0(wd00,"/data/MONIS6_2021_data")
wd00.1 <- wd00_wddata
wd00.1_inf04 <- paste0(wd00.1,"/",inf04 )
# Read in the xlsx file
df_dtc_asss <- openxlsx::read.xlsx(wd00.1_inf04,1)
# combine the columns with the species and the abbreviations
# to get columns with full species names and genus names and taxonomical
# higher levels
df_e13 <- dplyr::left_join(df_e12,
                df_dtc_asss, 
                by = c("speciesabbr" = "AbbrvNm"))

##########################################################################################
##########################################################################################
# section to try out different factors to adjust the resolution and to try out
# different factors to adjust the mean  neighbouring distance
# start --
##########################################################################################
#r_mnd <- seq(0.4,2, 0.1)
nr_mnd <- seq(0.2,2, 0.1)
res_m <- seq(0.01,0.5, 0.04)
res_m <- 0.01
res_m <- 0.06

nr_mnd <- 1.8
r_mnd <- 1.8
#mean neighbouring distance of 0.9 allows interpolation for
#Prorocentrum minimum
nr_mnd <- 0.9
r_mnd <- 0.9
#res_fac <- 0.25
#res_m <- 0.25
##########################################################################################
# section to try out different factors to adjust the resolution and to try out
# different factors to adjust the mean  neighbouring distance
# end --
##########################################################################################

##########################################################################################
# based on the section above trying out different factors to adjust
# the resolution and to try out the mean  neighbouring distance - use this
# data frame for the species
# start --
##########################################################################################
u_spc_yr <-   df_e13 %>% dplyr::distinct(Latinsk_navn, Dato_inds.yy)

#unq.spc.years <- "Mnemiopsis leidyi.2018"
# col.h2 <- c("spc.lat.nm","mean.neigh.dist","res_fact")
# Mnelei <- c("Mnemiopsis_leidyi", 1.2, 0.25) # coarse resolution = faster to calculate
# #Mnelei <- c("Mnemiopsis_leidyi", 1.2, 0.05) # fine resolution = slower to calculate
# Myaare <- c("Mya_arenaria",0.4, 0.25) # coarse resolution = faster to calculate
# #Myaare <- c("Mya_arenaria",0.4, 0.03) # fine resolution = slower to calculate
# Colper <- c("Colpomenia_peregrine",0.5,0.25) # coarse resolution = faster to calculate
# #Colper <- c("Colpomenia_peregrine",0.01,0.25) # fine resolution = slower to calculate
# Psever <- c("Pseudochattonella_verruculosa",1.1,0.25)
# Psefar <- c("Pseudochattonella_farcimen",0.5,0.25) # coarse resolution = faster to calculate
# Karmik <- c("Karenia_mikimotoi",0.35,0.25)
# Bonham <- c("Bonnemaisonia_hamifera", 0.4, 0.25)
# Promin <- c("Prorocentrum_minimum", 0.8,0.25)
# Cragig <- c("Crassostrea_gigas",0.4,0.1)
#
#
# col.h2 <- c("spc.lat.nm","res_fact","mean.neigh.dist")
# Acibae <- c("Acipenser_baerii",0.05,0.5)
# Bonham <- c("Bonnemaisonia_hamifera",0.05,0.5)
# Caraur <- c("Carassius_auratus",0.05,0.5)
# Colper <- c("Colpomenia_peregrine",0.01,0.5)
# Corcas <- c("Cordylophora_caspia",0.05,0.5)
# Cragig <- c("Crassostrea_gigas",0.05,0.5)
# Cypcar <- c("Cyprinus_carpio",0.05,0.5)
# Erisin <- c("Eriocheir_sinensis",0.05,0.5)
# Homame <- c("Homarus_americanus",0.05,0.5)
# Karmik <- c("Karenia_mikimotoi",0.05,0.5)
# Mnelei <- c("Mnemiopsis_leidyi",0.05,0.5)
# Myaare <- c("Mya_arenaria",0.01,0.5)
# Neomel <- c("Neogobius_melanostomus",0.05,0.5)
# Oncmyk <- c("Oncorhynchus_mykiss",0.05,0.5)
# Oncgor <- c("Oncorhyncus_gorbuscha",0.05,0.5)
# Parcam <- c("Paralithodes_camtschaticus",0.05,0.5)
# Promin <- c("Prorocentrum_minimum",0.09,0.5)
# Psefar <- c("Pseudochattonella_farcimen",0.05,0.5)
# Psever <- c("Pseudochattonella_verruculosa",0.05,0.5)
# Rhihar <- c("Rhithropanopeus_harrisii",0.05,0.5)
#
# #fc.spc_df <- as.data.frame(rbind(Mnelei, Myaare, Colper, Psever, Psefar, Karmik, Bonham, Promin, Cragig))
# fc.spc_df <- as.data.frame(rbind(Acibae, Bonham, Caraur, Colper, Corcas, Cragig, Cypcar, Erisin, Homame, Karmik, Mnelei, Myaare, Neomel, Oncmyk, Oncgor, Parcam, Promin, Psefar, Psever, Rhihar))
# colnames(fc.spc_df) <- col.h2
# #get first two rows of df
# #fc.spc_df <- fc.spc_df[1:2,]

##########################################################################################
# based on the section above trying out different factors to adjust
# the resolution and to try out the mean  neighbouring distance - use this
# data frame for the species
# end --
##########################################################################################


##########################################################################################
# use this loop below to try out different factors to adjust
# the resolution and to try out the mean  neighbouring distance
# start --
##########################################################################################

#for (r_mnd in nr_mnd) { #loop for mean neighbouring distance
for (res_fac in res_m) { #loop for resolution in costras
  print(paste("mnd:",r_mnd))
}

##########################################################################################
# use this loop below to try out different factors to adjust
# the resolution and to try out the mean  neighbouring distance
# end -- NOTICE that the curly bracket for this loop is ending further down !!
##########################################################################################

#################################################################################
# BEGIN LOOP over species for making 2 ipdw maps , one map per season
#################################################################################
#Extract Unique Elements from main data frame
#unq.spc.seas <- unique(loc_edna01$spc.season)
df_e13$Latinsk_navn <- gsub("^ ","",df_e13$Latinsk_navn)
df_e13$Latinsk_navn <- gsub(" $","",df_e13$Latinsk_navn)
df_e13$LatNm_wu <- gsub(" ","_",df_e13$Latinsk_navn)
# get the unique species names
unq.spc <- unique(df_e13$Latinsk_navn)
#replace underscore with space
unq.spc <- gsub("_"," ",unq.spc)

#uncomment below if you want to test the iteration
# with a single species
#unq.spc <- "Colpomenia_peregrine"
#unq.spc <- "Mnemiopsis leidyi"
#unq.spc <- "Pseudochattonella farcimen"
#unq.spc <- "Mya arenaria"
#unq.spc <- "Bonnemaisonia hamifera"
#unq.spc <- "Bonnemaisonia hamifera"
#unq.spc <- "Crassostrea gigas"
#unq.spc <- "Pseudochattonella verruculosa"
#unq.spc <- "Prorocentrum minimum"
#Extract Unique Elements from shortened dataframe
yrs <- unique(as.numeric(df_e13$Dato_inds.yy))
#yrs <- "2018"
#count the elements
no.y2 <- length(yrs)
# make sequence of numbers
no.e1 <- seq(1:no.y2)
#bind columns and make a data frame
no.e2 <- as.data.frame(cbind(no.e1, yrs))
#get maximum number of years
nppy <- max(no.e1)

# make a vector with variables that represents the machines
# used for analysis
mchns <- unique(df_e13$mch)
# make a number that represents the number of machines
nmchns <- length(mchns)

# get the last element in the vector
unq.spc_n <- tail(unq.spc, n=1) 
# or try with all species names in the list
 unq.spc_n <- unq.spc
# -loop over species
for (spec.lat in unq.spc_n) {
  print(spec.lat)
  #}
  # get an index number for the species
  idxNosp <- which(spec.lat==unq.spc_n)
  # use the index number for the species appendix number
  no.spc.app.plot <- as.character(idxNosp)
  #get the latin species name without underscore
  spec.lat.w_undersc <- paste(sub(' ', '_', spec.lat))
  spec.lat.no_undersc <- spec.lat
  #subset the dataframe based on variable value in column
  df_e14 <- df_e13[which(df_e13$LatNm_wu == spec.lat.w_undersc),]
  #get the number for the appendix plot number
  #get the latin species name with an underscore
  spec.lat.w_undersc <- sub(' ', '_', spec.lat)
  #use two functions together 
  genusl <- substr(spec.lat.w_undersc, 1, 1)
  spcl <- gsub("*.*_","",spec.lat.w_undersc)
  #and paste them together
  # to get an abbreviated species name
  short.spec.lat <- paste(genusl,"_",spcl,sep="")
  
  
  #pad with zeros to 2 characters for 
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  TwN <- ifelse(nchar(idxNosp)<2,stringr::str_pad(idxNosp, 2, pad = "0"),idxNosp)
  TwN <- as.character(TwN)
  sbs.AssIDNo <- TwN
  # use the data frame with best inferred resolutions and multiplier for
  # mean neighbouring distance
  # These factors were inferred using sequences of different multipliers
  # to see which returned a detailed heat map for each of the species
  #match species name to get resolution factor
  #res_fac <- fc.spc_df$res_fact[match(spec.lat, fc.spc_df$spc.lat.nm)]
  #res_fac<-as.numeric(as.character(res_fac))
  #match species name to get  factor for multiplying mean neighbouring distance
  #r_mnd <- fc.spc_df$mean.neigh.dist[match(spec.lat, fc.spc_df$spc.lat.nm)]
  #r_mnd<-as.numeric(as.character(r_mnd))
  # Exporting PFD files via postscript()
  #first make a file name
  flnm <- paste(
    wd00_wd13,"/Fig13_v",sbs.AssIDNo,
    "_ipwd_res",res_fac,"_mnd",
    r_mnd,
    "_",
    short.spec.lat,
    ".png",
    sep = ""
  )
  #use the filename in the plot to print to pdf
  png(file=flnm,
      #define the dimensions on the page in the pdf file
      width = ((2) *297),
      height = ((nppy*nmchns)* 210),
      units = "mm",
      # define the resolution of the plot , you must add resolution
      # otherwise the plot will not be produced
      res = 300)
  # if this throws an error about not being able to open the pdf it is most
  # likely because you have the Adobe Acrobat reader open already.
  # It cannot write the file to a file you already have open in another application
  # Close the Adobe Acrobat reader, and try again
  #dev.off()
  #set plotting margins
  op <-
    par(
      mfrow = c(nppy*nmchns, 2), # set number of panes inside the plot 
      # - i.e. c(2,2) would make four panes for plots 
      # - "c('plots in rows', 'plots in columns')"
      oma = c(1, 1, 0, 0), # set outer margin (the margin around the combined plot area) - higher numbers increase the number of lines
      mar = c(5, 5, 3, 5) # set the margin around each individual plot . with sides ordered as: " c(bottom,left,top,right)"
    )
  #dev.off()
  # to test the loop with one variable  use a single selected year
  #yrs <- "2017"
  #yr_smpl <- "2017"
  
  #loop over years sampled -  to produce individual tables per year sampled
  for (yr_smpl in yrs){
    print(yr_smpl)
    #}
    
    #get maximum eDNA log 10 level
    #get prefered upper limit of legend, rounded up - use later for zlim
    zvm<-ceiling(max(df_e14$l10cp_L))
    
    #subset the original dataframe per year sampled
    #to be able to plot eDNA evaluation squares
    df_e15 <-
      df_e14[which(df_e14$Dato_inds.yy == yr_smpl),]
    # to test the loop with one variable  use a single selected season
    # season <- "season_2"
    #categories.of.seasons <-"season_1"
    #reverse the order of the elements in the vector
    # because you will want to have the plot with the ipdw
    #map for spring on the left, and the fall ipdw map on the right
    rctse <- rev(categories.of.seasons)
    #iterate over the season in the vector
    for (season in rctse){
      print(season)
      #}
      
      #iterate over the machines in the vector
      for (mch_tp in mchns){
        print(mch_tp)
        #}
        
      #use match to match the season with a data frame and get the name for the season
      spcfc_seaon_name <- seaons_nms_df$names.of.seasons[match(season, 
                                    seaons_nms_df$categories.of.seasons)]
      spcfc_seaon_name <- as.character(spcfc_seaon_name)
      
      # subset in the  data frame to get the data for the season
      # and the machine type
      df_e15_ssn <- df_e15[ which((df_e15$ssn.smpl== season & df_e15$mch== mch_tp)),]
      
      
      #check if the data frame has more than 2 rows
      # See this question : https://stackoverflow.com/questions/35366187/how-to-write-if-else-statements-if-dataframe-is-empty
      #check across two statements
      # if the data frame has less then 2 rows, then it will not be possible
      # to to the interpolation, and then it needs to be skipped
      if (dim(df_e15_ssn)[1] <= 2) {
        print(paste("data frame for",spcfc_seaon_name,yr_smpl,"has less than 2 rows",
                    "and cannot be interpolated on",sep=" "))
        #}
        #if subsetted data fram for spring is empty - no samples
        #then create an empty data frame with zeroes and no color for points
        plot.new()
        #dev.off()
      }
      # check for multiple situations - second check if the data
      # frame has more than one dimension - 
      # https://www.datamentor.io/r-programming/if-else-statement/
      # if the dataframe does have more than one dimension, then try and plot it
      else if (dim(df_e15_ssn)[1] >= 3)
        #start curly bracket for 'else if' test testing 
        # whether the dimensions on the data frame is more than one, 
        # if it is then make the plot on the map
      {
        
        #get minimum and maximum to define range for jitter of points
        M27_jmin_lon <- min(df_e15_ssn$lon.m)
        M27_jmax_lon <- max(df_e15_ssn$lon.m)
        jit_lon <- (M27_jmax_lon-M27_jmin_lon)/80000
        #get minimum and maximum to define range for jitter of points
        M27_jmin_lat <- min(df_e15_ssn$lat.m)
        M27_jmax_lat <- max(df_e15_ssn$lat.m)
        jit_lat <- (M27_jmax_lat-M27_jmin_lat)/80000
        #jitter points to work around overlapping points
        df_e15_ssn$jit.lok_pos_lon <- jitter(df_e15_ssn$lon.m, jit_lon)
        df_e15_ssn$jit.lok_pos_lat <- jitter(df_e15_ssn$lat.m, jit_lat)
        #make SpatialPointsDataFrame with decimal degree coordinates
        #make points with decimal coordinates - this makes use of a different CRS
        pnts2 <- SpatialPointsDataFrame(df_e15_ssn[,c("jit.lok_pos_lon",
                                                      "jit.lok_pos_lat")],
                                        df_e15_ssn,
                                        proj4string = crs(epsg4326nCRS2))
        #make another data frame
        pnts3 <- SpatialPointsDataFrame(df_e15_ssn[,c("lon.m",
                                                      "lat.m")],
                                        df_e15_ssn,
                                        proj4string = crs(epsg4326nCRS2))
        #pnts_test <- as.data.frame(pnts2)
        #make points in decimal degr coordinates in a spatial dataframe for sampl locations
        #but moved a bit aside to plot eDNA evaluations - for the season
        df_e15_ssn$lok_pos_lon.f.pnts4 <- df_e15_ssn$jit.lok_pos_lon + 0
        df_e15_ssn$lok_pos_lat.f.pnts4 <- df_e15_ssn$jit.lok_pos_lat + 0.16
        
        
        #turn these new points in to a spatial data frame
        pnts4 <-
          SpatialPointsDataFrame(df_e15_ssn[, 
                                           c("lok_pos_lon.f.pnts4", 
                                             "lok_pos_lat.f.pnts4")],
                                 df_e15_ssn,
                                 proj4string = crs(epsg4326nCRS2))
        
        #prepare bounding box for decimal degrees
        bbox_k2 <- raster::buffer(
          as(extent(sp::spTransform(pnts2, projection(coastline10))), 
             "SpatialPolygons"),
          width = 7) # A width=10 zooms up and includes more 
        # sourrounding landmass  in the map. A width=2 zooms in and 
        # includes less landmass around in the map.
        #project coastline on bbbox
        projection(bbox_k2) <- projection(coastline10)
        
        # The 'raster::crop' does not work any longer
        #pols2 <- raster::crop(coastline10, bbox_k2)
        # following the advice here instead
        #https://stackoverflow.com/questions/78242946/error-unable-to-find-an-inherited-method-for-function-crop-for-signature-s
        library(sf)
        st_is_valid(coastline10, reason = TRUE)[!st_is_valid(coastline10)]
        #transform to spdf with decim degr coordinate
        sf_use_s2(use_s2 = FALSE)
        #> Spherical geometry (s2) switched off
        csl10_crop <- st_crop(coastline10, bbox_k2)
        csl10_crop <- st_crop(coastline10, xmin=6, xmax=17, ymin=54, ymax=59)
        pols2 <- csl10_crop
        #> although coordinates are longitude/latitude, st_intersection assumes that they
        #> are planar
        #> Warning: attribute variables are assumed to be spatially constant throughout
        #> all geometries
        sf_use_s2(use_s2 = TRUE)
        #> Spherical geometry (s2) switched on
        # plot(csl10_crop$geometry)
        # dev.off()
        #transform to spdf with decim degr coordinates
        # if 'pols2' is a data frame, then transform it to a spatial object
        # using 'st_transform' otherwise use 'sp::spTransform'
        # https://gis.stackexchange.com/questions/307432/sptransform-error-unable-to-find-inherited-method-for-function-sptransform
        if (is.data.frame( pols2)==T){
          pols2 <- st_transform(pols2, projection(pnts2))  
        } else {
          pols2 <- sp::spTransform(pols2, projection(pnts2))
        }
        #plot(pols2)
        
        #make costras with lon lat  coordinates
        costras2 <- costrasterGen(pnts2, pols2, extent = "polys",
                                  projstr = projection(pols2),
                                  resolution = res_fac)
        
        #resolution = 0.08) #Note this is in decimals for lonlat 
        #coordinates. Higher value (>1) makes interpolation goes fast, 
        #low (<0.05) takes too much RAM and too much time
        #make ipdw result 'res.ipdw2' based on pnts with lon-lat 
        #coordinates
        #This ipdw result 'res.ipdw2' is different in blending between 
        #colours between points
        #As compared to the section below that involves the 'training' 
        #element in 'res.ipdw3'
        # res.ipdw2 <- ipdw::ipdw(pnts2, 
        # costras2, paramlist = "copy.per.L.log10",
        #                         range = 10, 
        # low range reduces areas around point. 
        # Note this is in decimals for lonlat coordinates.
        #                         dist_power = 4)
        # increasing the 'dist_power' to 10 makes gradients between points 
        # less pronounced, and instead colors more intensely up to the borders
        
        ###################################################################################
        pnts <- pnts2
        costras <- costras2
        # find average nearest neighbor
        library(spatstat)
        
        W              <- owin(range(coordinates(pnts)[,1]),
                               range(coordinates(pnts)[,2]))
        kat.pp         <- ppp(coordinates(pnts)[,1], 
                              coordinates(pnts)[,2], window = W)
        #if the ppp function complains about duplicated points, you can check which ones are duplicated
        #duplicated(kat.pp)
        #go back and check your original lon-lat columns in the input data frame and make sure no lon-lat are duplicated.
        mean.neighdist <- mean(nndist(kat.pp))
        # grid building
        gridsize       <- mean.neighdist * 1*r_mnd #increasing the multiplier makes the gradients blend more across the borders equidistant from each point
        r_mnd1.lonlat       <- (mean.neighdist * 1*r_mnd)
        grainscale.fac <- gridsize / res(costras)[1]
        gridras        <- aggregate(costras, fact = grainscale.fac)
        #as.data.frame(costras)
        gridpol        <- rasterToPolygons(gridras)
        gridpol$value  <- row.names(gridpol)
        #check pnts as a dataframe
        #df.pnts01  <- data.frame(pnts)
        #colnames(df.pnts01)
        # spatial join
        fulldataset.over    <- over(pnts, gridpol)
        fulldataset.over    <- cbind(data.frame(fulldataset.over),
                                     setNames(data.frame(pnts),
                                              c(colnames(data.frame(pnts)))))
        #colnames(data.frame(pnts))
        #fulldataset.over
        # grid selection
        set.seed(2)
        gridlev <- unique(fulldataset.over$value)
        for(i in seq_along(gridlev)){
          activesub <- subset(fulldataset.over, fulldataset.over$value == gridlev[i])
          selectnum <- gdata::resample(seq_len(nrow(activesub)), 1)
          if(i == 1){
            training <- activesub[selectnum,]
          }
          else{
            training <- rbind(training, activesub[selectnum,])
          }
        }
        #####
        validate             <- fulldataset.over[!(row.names(fulldataset.over) %in%
                                                     row.names(training)),]
        #Make sure you change the 'lon_decp02' and 'lat_decp02' to your own lon-lat columns in your original dataframe
        #Make sure you change the 'jit.lok_pos_lon' and 'jit.lok_pos_lat' to your own lon-lat columns in your original dataframe
        xy                   <- cbind(training$jit.lok_pos_lon, training$jit.lok_pos_lat)
        training             <- SpatialPointsDataFrame(xy, training)
        xy                   <- cbind(validate$jit.lok_pos_lon, validate$jit.lok_pos_lat)
        validate             <- SpatialPointsDataFrame(xy, validate)
        projection(training) <- projection(pnts)
        projection(validate) <- projection(pnts)
        ####
        
        #paramlist is the z-value to interpolate across
        # use "log.10_copies_L" as input parameter
        # to use eDNA lvls that at least are above zero 
        pl <- c("log.10_copies_L")
        # use the "cpLwp1l10" to only use eDNA lvls 
        # above LOQ
        pl <- c("cpLwp1l10") # use the copy number per L plus 1
        pl <- c("l10cp_L") # use the copy number per L plus 1
        #but only for above LOQ
        class(training)
        # convert to sf object
        training_sf <- st_as_sf(training)
        class(training_sf)
        #make the ipdw raster
        res.ipdw3 <- ipdw::ipdw(training_sf, costras, range = mean.neighdist * 8*r_mnd, pl,
                          overlapped = TRUE, dist_power = 1.0)
        r_mnd2.lonlat <- (mean.neighdist * 8 *r_mnd)
        #low range reduces areas around point. #Note this is in decimals for lonlat coordinates.
        #increasing the 'dist_power' to 10 makes gradients between points less pronounced, and instead colors more intensely up to the borders
        
        
        #get color from the subsetted data frame
        col.f.ramp.pal <- "blue" # unique(loc_edna04$col_f_Phyl)
        col.f.ramp.pal <- "orange" #unique(loc_edna04$spcf_col_f_Phyl)
        #match species name in loop with the color for the species
        #phy_col <- tx_hierc_df$col_f_Phyl[match(spec.lat,tx_hierc_df$genus_spec)]
        phy_col <- "darkgreen"
        #okace the color in the object that is needed in the plot
        col.f.ramp.pal <- phy_col
        #col.f.ramp.pal <- "darkgreen" #unique(loc_edna04$spcf_col_f_Phyl)
        #make a colour range
        colfunc03 <-
          colorRampPalette(c("white", col.f.ramp.pal, "black"))
        Apts = 30
        #make color range
        cols01 <- colfunc03(30) #white to phylum-color to black
        #get prefered upper limit of legend, rounded up - use later for zlim
        #zvm<-ceiling(max(loc_edna04$copy.per.L.log10))
        #define limits to the axis of the plot
        xx <- labeling::extended(3, 17, 14, only.loose=TRUE)
        yy <- labeling::extended(54, 60,6, only.loose=TRUE)
        zz <- labeling::extended(0, zvm,zvm, only.loose=TRUE)
        #define limits to use to crop map
        e <- extent(3, 17, 54, 59.5)
        #crop the map before plotting
        res.ipdw3c <- crop(res.ipdw3, e)
        
        #plot heatmap 02 - with decimal degress and lon lat
        plot(res.ipdw3c, #main = "copy.per.L.log10 lonlat train",
             col=cols01, #set color for z-value
             zlim=c(0,zvm), #limit the z-value in the legend and map
             xlim=c(3,17),
             ylim=c(54,60),
             xaxt="n", #do not use the default tick labels on the x-axis
             yaxt="n", #do not use the default tick labels on the y-axis
             #xlab = "Longitude", #label along x-axis
             #ylab = "Latitude", #label along y-axis
             las=1,  #las=1 turns the tick label on the axis
             cex.axis=2) #adjust size
        #Add text along x- and y-axis in specified color and size
        mtext("longitude", side=1, line=3, col="black", cex=2)
        mtext("latitude", side=2, line=3, col="black", cex=2)
        #add land
        plot(pols2, add = TRUE, col="black", bg="azure4")
        #add land in a different colour
        plot(pols2, add = TRUE, col="azure4")
        
        #add the the eDNA evaluations for sampling points as symbols coloured by column
        plot(pnts4,add = TRUE,pch = 22,col = "blue",
             bg = c(as.data.frame(pnts4)$eDNA_eval_t_repl_col),cex = 3.0)
        # add sampling locations coloured by conventional monitoring result
        # plot(pnts3,add = TRUE,pch = 25,col = "blue",bg = c(as.data.frame(pnts3)$col.f.conv_rec_val),cex = 2.6)
        #add special ticks and adjust the size of the text associated 
        #with these tick marks
        axis(1, at = xx, cex.axis=2.0)
        axis(2, at = yy, las=1, cex.axis=2.0) #las=1 turns the labels at the tick marks on the axis
        #set range for z key
        zz <- seq(0,zvm)
        #axis(5, at = zz, las=1, cex.axis=2.0) #las=1 turns the labels at the tick marks on the axis
        #deduct a bit from the latitude, to lower the positioning of 
        #the label
        #redlatpt1 <- (max(as.data.frame(pnts4)$declat)-min(as.data.frame(pnts4)$declat))/350.6
        #add text to the same lon lat position for harbour, or for month
        
        # text(as.data.frame(pnts4)$declon,
        #      as.data.frame(pnts4)$declat-redlatpt1,
        #      as.data.frame(pnts4)$Harbour.abbr1,
        #      #as.data.frame(pnts4)$Harbour,
        #      #as.data.frame(pnts4)$month2, # use this one instead if you want the sampling month underneath
        #      pos=4,
        #      cex= 1.8
        # )
        # 
        #add a title for key-legend on side
        mtext("log10(eDNA copies/L)", side=4, line=5, cex=1.8)
        latnm <- spec.lat.w_undersc
        spec.lat.no_undersc <- paste(sub('_', ' ', latnm))
        #add a title with bquote
        #use 'atop' function to get the title to appear on two lines
        #see this website: https://stackoverflow.com/questions/20549337/expression-and-new-line-in-plot-labels?lq=1
        title(main = c(bquote(
          atop(
            'eDNA lvls-log10 from'
            ~ italic(.(spec.lat.no_undersc)),
            ~ '('
            ~ .(mch_tp)
            ~ '), Assay Id No:' ~ .(sbs.AssIDNo)
            ~ ' for '
            ~ .(spcfc_seaon_name)
            ~ ' in '
            ~ .(yr_smpl)
            ~ '  ')
        )))
        # ##add legend for conventional monitoring  evaluation
        # legend(
        #   "bottomleft",
        #   "(x,y)",
        #   bg = "white",
        #   c(
        #     "Unknown, no previous record",
        #     "Recorded in the past",
        #     "Recorded in 2017"
        #   ),
        #   ncol = 1,
        #   pch = c(25, 25, 25),
        #   pt.cex=2.2, #size of points in legend
        #   #set type of point
        #   #col= c("black", "black", "black"), #set color of point
        #   col = c("blue", "blue", "blue"),
        #   #set color of point
        #   #pt.bg=c(alpha(c("white", "red", "black"), 0.6)), #set background color of point
        #   pt.bg = c(c("white", "red", "black")),
        #   #set background color of point
        #   pt.lwd = c(1.0),
        #   title = "conventional monitoring",
        #   cex = 1.4,
        #   inset = 0.02
        # )
        
        # # add legend for spring and fall
        # legend(
        #   "topleft",
        #   "(x,y)",
        #   bg = "white",
        #   c("spring", "autumn"),
        #   ncol = 1,
        #   pch = c(22, 22),
        #   #set type of point
        #   col = c("blue", "red"),
        #   #set color of point
        #   #pt.bg=c(alpha(c("white", "white"), 0.6)), #set background color of point
        #   pt.bg = c(c("white", "white")),
        #   #set background color of point
        #   pt.lwd = c(1.2),
        #   title = "season ",
        #   cex = 1.1,
        #   inset = 0.02
        # )
        
        # # # add legend for eDNA evaluation
        # legend(
        #   "topright",
        #   "(x,y)",
        #   bg = "white",
        #   c(
        #     "No Ct",
        #     "below LOD",
        #     "above LOD and below LOQ" ,
        #     "1 above LOQ",
        #     "3 above LOQ"
        #   ),
        #   ncol = 1,
        #   pch = c(22, 22, 22, 22, 22),
        #   pt.cex=2.4,
        #   #set type of point
        #   col = c("black", "black", "black", "black", "black"),
        #   #set color of point
        #   #pt.bg=c(alpha(c("white", "yellow", "black"), 0.6)), #set background color of point
        #   pt.bg = c(c(
        #     "white", "yellow", "orange", "red", "black"
        #   )),
        #   #set background color of point
        #   pt.lwd = c(1.0),
        #   title = "eDNA evalution ",
        #   cex = 1.4,
        #   inset = 0.02
        # )
        
        
        #end curly bracket for 'else if' test testing whether the dimensions on the data frame is more than one, if it is then make the plot on the map
      }
      
      #end loop over machine
      }
      
      #end loop over seasons
    }
    #end loop over year
  }
  
  
  
  # add title for the pdf-page
  mtext(
    c(paste("Fig13_v",idxNosp,"_Appendix E", no.spc.app.plot, "."),  sep = ""),
    outer = TRUE,
    #use at , adj and padj to adjust the positioning
    at = (par("usr")[1] + 0.15 * diff(par("usr")[1:2])),
    adj = 3.4,
    padj = 2,
    #use side to place it in the top
    side = 3,
    cex = 1.6,
    line = -1.15)
  #apply the par settings for the plot as defined above.
  par(op)
  # end pdf file to save as
  dev.off()
  
  ##end loop over species
}
##

#################################################################################
# END LOOP over species for making 2 ipdw maps , one map per season
#################################################################################

# use this end of loop to loop over factors adjusting the resolution in the costras
# and use the same loop to try out different factors adjusting the mean neighbouring distance
#}






#################################################################################


####################################################################################
# End Appendix E
####################################################################################


dev.off()

#define output file and path
outfile01 <- paste(wd00,wd13,"/suppmatr_13_01_MONIS4eDNA01_df.csv",sep="")
#write the table to a csv
#to be used for the next R-code that plots eDNA levels on maps
write.csv(MONIS4eDNA01_df, file = outfile01)

#MONIS4eDNA01_df

