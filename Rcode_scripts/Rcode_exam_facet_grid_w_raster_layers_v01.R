#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

library(dplyr)
library(tidyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspatial)
library(patchwork)
library(cowplot)
#dev.off()
#____________________________________________________________________________#
# R-code provided for the project:
# “comparison of qPCR and ddPCR”
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
# #Delete any previous versions of the output directory
# unlink(wd00_wd13, recursive=TRUE)
# #Create a directory to put resulting output files in
# dir.create(wd00_wd13)


wd06 <- "output06_compare_ddPCR_and_qPCR_results"
wd00_wd06 <- paste(wd00,"/",wd06,sep="")
inf05 <- paste0(wd00_wd06,"/","df_csvs05.2_v01.csv")
df_csvs05 <- read.table(inf05,header=T, sep = ",")
#______________________________
# define directory and read in the table w colors for species
wdd <- "data/MONIS6_2021_data"
wd00_wdd <- paste0(wd00,"/",wdd)
fl.in <- "table_w_colors_for_spc.csv"
fl.in <- paste0(wd00_wdd,"/",fl.in)
df_clr2 <- read.table(fl.in,header=T, sep = ";")

#################################################################################

# define folder with previously stored csv file
wdout11 <- "output11_compare_all_ddpcr_and_qpcr"
flin <- paste0(wd00,"/",wdout11,
               "/table_w_ddpcr_and_qpcr_data_and_lat_lon_v03.csv")

# read in the data frame as a csv file
df_e12 <- read.csv(file=flin,
                   check.names = F,
                   fileEncoding = "Latin1",
                   header =T)
# get the column headers
clhd12 <- colnames(df_e12)
# replace to only ASCII characters in the column headers
clhd12 <- iconv(clhd12, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
# also use the stringi package to replace to ASCII characters
clhd12 <- stringi::stri_trans_general(clhd12, "latin-ascii")
clhd12 <- iconv(clhd12, "latin1", "ASCII", sub="")
clhd12 <- gsub("[^\u0001-\u007F]+|<U\\+\\w+>","", clhd12)
clhd12 <- textclean::replace_non_ascii(clhd12)
colnames(df_e12) <- clhd12
# identify the column headers that ends with '.y'
cl12.tom<- which(grepl("\\.y$",clhd12))
# omit these columns from the data frame
df_e12 <- df_e12[,-cl12.tom]
# get all the column names and order them alphabetically
clhd12 <- colnames(df_e12)
clhd12 <- clhd12[order(clhd12)]
# define a vector with columns to exclude from the data frame
cl.to <- c("U_Pr_Nr.x","Lok_omr01.x","lokalitet_vanda.x" , "Dato_inds.x")
df_e12 <- df_e12[,!(names(df_e12) %in% cl.to)]
# use gsub to remove '.x' the the ending of all column headers
colnames(df_e12) <- gsub("\\.x$","", colnames(df_e12))


# use left join to add the colors for each species
df_e12 <- dplyr::left_join(df_e12,
                           df_clr2, 
                           by = c("speciesabbr" = "spcNm"))
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
# subset the 'df_dtc_asss' data frame  to only have distinct
# rows for 'Genus' and 'species, and keep all colunms
df_dtc_asss <- dplyr::distinct(df_dtc_asss, Genus, Species, .keep_all = TRUE)
# use gsub to replace leading white space in the column 'Latinsk_navn'
df_dtc_asss$Latinsk_navn <- gsub("^ ","",df_dtc_asss$Latinsk_navn)
# exclude the column 'Latinsk_navn' from the 'df_dtc_asss' data frame
cltom <- which(grepl("Latinsk_navn",colnames(df_dtc_asss)))
df_dtc_asss <- df_dtc_asss[,-cltom]
#df_dtc_asss$AbbrvNm
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

df_e13$eDNA_eval_t_repl_col <- "red"

colnames(df_e13)
#View(df_e13)

df_e13$speciesabbr
df_e13$mch
df_e13$smplNm
df_csvs05 %>% dplyr::select(spcAbr)
#View(df_csvs05)
df_csvs05$mch <- df_csvs05$vlc
df_csvs05$speciesabbr <- df_csvs05$spcAbr
df_csvs05$smplNm
#identify colnames identical in both df_csvs05 and df_e13
cl_idn <- intersect(colnames(df_csvs05), colnames(df_e13))
# exclude the "eDN.lv1.l10" column from the df_csvs05
df_csvs06 <- df_csvs05 %>% dplyr::select(-eDN.lv1.l10)
# use left_join to join the two data frames
# using the columns 'speciesabbr' , 'mch' and 'smplNm'
# as common columns
df_e13 <- dplyr::left_join(df_e13,
                           df_csvs06,
                           by = c("speciesabbr","mch","smplNm"))
#View(df_e13)
# make a vector with variables that represents the machines
# used for analysis
mchns <- unique(df_e13$mch)
# make a number that represents the number of machines
nmchns <- length(mchns)
#dev.off()
# get the last element in the vector
unq.spc_n <- tail(unq.spc, n=1) 
#dev.off()
#unq.spc_n <- unq.spc_n[16]
# #
sp.t.incl <- which(grepl("verruculosa",unq.spc_z))
unq.spc_n <- unq.spc_z[sp.t.incl]
# str(unq.spc_n)
# str(df_e13)

org.df_e13 <- df_e13
df_e13 <-org.df_e13

df_e13 <- df_e13 %>%
  dplyr::select(col.f_spc,
                Latinsk_navn,
                LatNm_wu,
                l10cp_L,
                Dato_inds.yy,
                ssn.smpl,
                mch,
                lon.m,
                lat.m,
                cl.f.evl)
# subset the dataframe to only include rows
# where the species names contains "Pseudochattonella verruculosa"
#unq.spc_n <- unique(df_e13$Latinsk_navn)
#df_e13 <- df_e13[which(df_e13$Latinsk_navn %in% "Colpomenia peregrine"), ]
#df_e13 <- df_e13[which(df_e13$Latinsk_navn %in% "Pseudochattonella verruculosa"), ]
#str(df_e13)
#View(df_e13)
# -loop over species
#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspatial)
library(patchwork)
library(cowplot)
library(ipdw)
library(scales)
library(sf)
library(rnaturalearth)
library(spatstat)
library(sp)

# [Previous code remains unchanged until the loop over species]

# BEGIN LOOP over species for making 2 ipdw maps , one map per season
# [Previous setup code remains unchanged]


# or try with all species names in the list
unq.spc_n <- unq.spc
sp.t.excl <- which(grepl("mykiss",unq.spc_n))
# exclude this species
unq.spc_n <- unq.spc_n[-sp.t.excl]
unq.spc_z <- unq.spc_n

# sp.t.incl <- which(grepl("peregrine",unq.spc_z))
# unq.spc_n <- unq.spc_z[sp.t.incl]

for (spec.lat in unq.spc_n) {
  print(spec.lat)
  # [Previous code for setting up spec.lat, clf.spc, etc.]
  print(spec.lat)
  # get the color for the species
  clf.spc <- unique(df_e13$col.f_spc[match(spec.lat, df_e13$Latinsk_navn)])
  # get an index number for the species
  idxNosp <- which(spec.lat == unq.spc_n)
  # use the index number for the species appendix number
  no.spc.app.plot <- as.character(idxNosp)
  # pad with zeros to 2 characters for 'no.spc.app.plot'
  no.spc.app.plot <- ifelse(nchar(no.spc.app.plot) < 2, stringr::str_pad(no.spc.app.plot, 2, pad = "0"), idxNosp)
  no.spc.app.plot <- as.character(no.spc.app.plot)
  # get the latin species name without underscore
  spec.lat.w_undersc <- paste(sub(' ', '_', spec.lat))
  spec.lat.no_undersc <- spec.lat
  spc_label <- spec.lat
  # subset the dataframe based on variable value in column
  df_e14 <- df_e13[which(df_e13$LatNm_wu == spec.lat.w_undersc), ]
  # use two functions together
  genusl <- substr(spec.lat.w_undersc, 1, 1)
  spcl <- gsub(".*_", "", spec.lat.w_undersc)
  # and paste them together
  short.spec.lat <- paste(genusl, "_", spcl, sep = "")
  # pad with zeros to 2 characters for
  TwN <- ifelse(nchar(idxNosp) < 2, stringr::str_pad(idxNosp, 2, pad = "0"), idxNosp)
  TwN <- as.character(TwN)
  sbs.AssIDNo <- TwN
  
  # Initialize a list to hold all plots for this species
  all_plots <- list()
  plot_counter <- 1
  # Initialize a list to hold all raster layers for this species
  
  lst_all_rst <- list()
  rast_cnt <-1
  # --- Calculate global max value for color scale ---
  zvm_global <- ceiling(max(df_e14$l10cp_L, na.rm = TRUE))
  # Initialize a list to hold all raster layers for this species
  lst_all_rst <- list()
  rast_cnt <- 1
  # --- Calculate global max value for color scale ---
  zvm_global <- ceiling(max(df_e14$l10cp_L, na.rm = TRUE))
  
  # loop over years sampled
  for (yr_smpl in yrs) {
    print(yr_smpl)
    # subset the original dataframe per year sampled
    df_e15 <- df_e14[which(df_e14$Dato_inds.yy == yr_smpl), ]
    # iterate over the season in the vector
    for (season in categories.of.seasons) {
      print(season)
      # [Previous code for setting up season, etc.]
      print(season)
      # use match to match the season with a data frame and get the name for the season
      spcfc_seaon_name <- seaons_nms_df$names.of.seasons[match(season, seaons_nms_df$categories.of.seasons)]
      spcfc_seaon_name <- as.character(spcfc_seaon_name)
      
      # iterate over the machines in the vector
      for (mch_tp in mchns) {
        print(mch_tp)
        # subset in the data frame to get the data for the season and the machine type
        df_e15_ssn <- df_e15[which((df_e15$ssn.smpl == season & df_e15$mch == mch_tp)), ]
        
        # Check if there are enough rows for interpolation
        if (nrow(df_e15_ssn) >= 3) {
          # --- ipdw interpolation ---
          # get minimum and maximum to define range for jitter of points
          M27_jmin_lon <- min(df_e15_ssn$lon.m)
          M27_jmax_lon <- max(df_e15_ssn$lon.m)
          jit_lon <- (M27_jmax_lon - M27_jmin_lon) / 80000
          # get minimum and maximum to define range for jitter of points
          M27_jmin_lat <- min(df_e15_ssn$lat.m)
          M27_jmax_lat <- max(df_e15_ssn$lat.m)
          jit_lat <- (M27_jmax_lat - M27_jmin_lat) / 80000
          # jitter points to work around overlapping points
          df_e15_ssn$jit.lok_pos_lon <- jitter(df_e15_ssn$lon.m, jit_lon)
          df_e15_ssn$jit.lok_pos_lat <- jitter(df_e15_ssn$lat.m, jit_lat)
          # make SpatialPointsDataFrame with decimal degree coordinates
          pnts2 <- SpatialPointsDataFrame(df_e15_ssn[, c("jit.lok_pos_lon", "jit.lok_pos_lat")], df_e15_ssn, proj4string = crs(epsg4326nCRS2))
          pnts3 <- SpatialPointsDataFrame(df_e15_ssn[, c("lon.m", "lat.m")], df_e15_ssn, proj4string = crs(epsg4326nCRS2))
          df_e15_ssn$lok_pos_lon.f.pnts4 <- df_e15_ssn$jit.lok_pos_lon + 0
          df_e15_ssn$lok_pos_lat.f.pnts4 <- df_e15_ssn$jit.lok_pos_lat + 0.16
          pnts4 <- SpatialPointsDataFrame(df_e15_ssn[, c("lok_pos_lon.f.pnts4", "lok_pos_lat.f.pnts4")], df_e15_ssn, proj4string = crs(epsg4326nCRS2))
          
          # --- Use coastline as barrier ---
          bbox_k2 <- raster::buffer(as(extent(sp::spTransform(pnts2, projection(coastline10))), "SpatialPolygons"), width = 7)
          projection(bbox_k2) <- projection(coastline10)
          library(sf)
          st_is_valid(coastline10, reason = TRUE)[!st_is_valid(coastline10)]
          sf_use_s2(use_s2 = FALSE)
          csl10_crop <- st_crop(coastline10, bbox_k2)
          csl10_crop <- st_crop(coastline10, xmin = 6, xmax = 17, ymin = 54, ymax = 59)
          pols2 <- csl10_crop
          sf_use_s2(use_s2 = TRUE)
          if (is.data.frame(pols2) == T) {
            pols2 <- st_transform(pols2, projection(pnts2))
          } else {
            pols2 <- sp::spTransform(pols2, projection(pnts2))
          }
          costras2 <- costrasterGen(pnts2, pols2, extent = "polys", projstr = projection(pols2), resolution = res_fac)
          pnts <- pnts2
          costras <- costras2
          library(spatstat)
          W <- owin(range(coordinates(pnts)[, 1]), range(coordinates(pnts)[, 2]))
          kat.pp <- ppp(coordinates(pnts)[, 1], coordinates(pnts)[, 2], window = W)
          mean.neighdist <- mean(nndist(kat.pp))
          gridsize <- mean.neighdist * 1 * r_mnd
          r_mnd1.lonlat <- (mean.neighdist * 1 * r_mnd)
          grainscale.fac <- gridsize / res(costras)[1]
          gridras <- aggregate(costras, fact = grainscale.fac)
          gridpol <- rasterToPolygons(gridras)
          gridpol$value <- row.names(gridpol)
          fulldataset.over <- over(pnts, gridpol)
          fulldataset.over <- cbind(data.frame(fulldataset.over), setNames(data.frame(pnts), c(colnames(data.frame(pnts)))))
          set.seed(2)
          gridlev <- unique(fulldataset.over$value)
          for (i in seq_along(gridlev)) {
            activesub <- subset(fulldataset.over, fulldataset.over$value == gridlev[i])
            # Skip if no rows in this grid cell
            if (nrow(activesub) == 0) next
            selectnum <- gdata::resample(seq_len(nrow(activesub)), 1)
            if (i == 1) {
              training <- activesub[selectnum, ]
            } else {
              training <- rbind(training, activesub[selectnum, ])
            }
          }
          validate <- fulldataset.over[!(row.names(fulldataset.over) %in% row.names(training)), ]
          xy <- cbind(training$jit.lok_pos_lon, training$jit.lok_pos_lat)
          training <- SpatialPointsDataFrame(xy, training)
          xy <- cbind(validate$jit.lok_pos_lon, validate$jit.lok_pos_lat)
          validate <- SpatialPointsDataFrame(xy, validate)
          projection(training) <- projection(pnts)
          projection(validate) <- projection(pnts)
          training_sf <- st_as_sf(training)
          pl <- c("l10cp_L")
          
          # Try to perform interpolation
          res.ipdw3 <- tryCatch({
            ipdw::ipdw(training_sf, costras, 
                       range = mean.neighdist * 8 * r_mnd, pl, 
                       overlapped = TRUE, dist_power = 1.0)
          }, error = function(e) {
            message("Interpolation failed for ",
                    spec.lat, " in ", season, 
                    " with ", mch_tp, ": ", e$message)
            NULL
          })
          
          # --- Plot with ipdw raster ---
          if (!is.null(res.ipdw3)) {
            # Convert raster to data frame for ggplot
            df_raster <- as.data.frame(res.ipdw3, xy = TRUE)
            names(df_raster) <- c("x", "y", "value") } else {
              #If the 'res.ipdw3' is NULL, then it is not possible
              # to adda a raster layer to the plot.
              # instead make a dummy data frame with NA values
              df_raster <- data.frame(x = numeric(0), 
                                      y = numeric(0),
                                      value = numeric(0))
            }
        } else {
          # If there are not enough points for interpolation
          # make a dummy data frame with NA values
          df_raster <- data.frame(x = numeric(0), 
                                  y = numeric(0),
                                  value = numeric(0))
        }
        # collect the raster layers into a list, with spec.lat, and year, and season,
        # and machine type in the name of each list element
        lst_all_rst[[rast_cnt]] <- list(
          df_raster = df_raster,
          season = season,
          mch_tp = mch_tp,
          yr_smpl = yr_smpl,
          spec.lat = spec.lat,
          zvm_global = zvm_global,
          clf.spc = clf.spc
        )
        # add to the raster counter
        rast_cnt <- rast_cnt + 1
      }
    }
  }
  #}
  # Now, prepare the data for plotting with the correct raster for each facet
  # We need to merge the raster data with the main data frame for faceting
  # So, we create a data frame that has all the info for each facet
  # and the corresponding raster data
  
  # Create a data frame for the facets
  facet_df <- expand.grid(
    ssn.smpl = categories.of.seasons,
    mch = mchns,
    stringsAsFactors = FALSE
  )
  
  
  # For each facet, find the matching raster in lst_all_rst
  facet_df <- facet_df %>%
    mutate(
      raster_data = purrr::map2(
        ssn.smpl, mch,
        ~ {
          match_idx <- sapply(lst_all_rst, function(x) {
            x$season == .x && x$mch_tp == .y
          })
          if (any(match_idx)) {
            lst_all_rst[[which(match_idx)]]$df_raster
          } else {
            data.frame(x = numeric(0), y = numeric(0), value = numeric(0))
          }
        }
      )
    )
  library(tidyr)
  library(tidyr)
  library(dplyr)
  
  # Ensure all combinations are present
  facet_df_expanded <- facet_df %>%
    expand(nesting(ssn.smpl, mch)) %>%
    left_join(facet_df, by = c("ssn.smpl", "mch"))
  
  
  # Replace missing raster_data with an empty data frame
  facet_df_expanded$raster_data <- lapply(facet_df_expanded$raster_data, function(x) {
    if (is.null(x)) {
      data.frame(x = numeric(0), y = numeric(0), value = numeric(0))
    } else {
      x
    }
  })
  
  # Unnest the raster data
  facet_df_unnested <- facet_df %>%
    unnest(cols = raster_data, keep_empty = T)
  #str(facet_df_unnested)
  colnames(facet_df_unnested)
  unique(facet_df_unnested$mch)
  
  
  pp <- ggplot() +
    theme_bw() +
    geom_sf(data = coastline10, fill = "white", color = "black") +
    geom_point(data = df_e14, aes(x = lon.m, y = lat.m, color = mch), size = 3, shape = 22) +
    geom_tile(
      data = facet_df_unnested,
      aes(x = x, y = y, fill = value),
      alpha = 0.7
    ) +
    scale_fill_gradientn(
      colors = c("white", clf.spc, "black"),
      values = scales::rescale(c(0, zvm_global / 2, zvm_global)),
      name = "log10(eDNA\ncopies/L)",
      limits = c(0, zvm_global)
    ) +
    coord_sf(xlim = c(6, 17), ylim = c(54, 59), expand = FALSE) +
    facet_grid(
      rows = vars(ssn.smpl),
      cols = vars(mch),
      drop = FALSE  # This ensures all facets are shown
    ) +
    theme(
      strip.background = element_rect(fill = "white", 
                                      colour = "white"),  # Change strip background and border color
      strip.text = element_text(face = "bold", hjust = 0, size=11),  # Make strip text bold and align to the left
      plot.title = element_text(hjust = -0.05, size = 16)  # Set title alignment and size
    ) +
    labs(
      title = bquote(bold("A. ") * italic(.(spc_label))),    # Mixed formatting for title
      x = "Longitude",  # Custom x-axis label
      y = "Latitude"    # Custom y-axis label
    )
  scale_color_manual(values = c("ddPCR" = "blue", "qPCR" = "red"))
  pp
  # Save the facet_grid plot
  flnm <- paste(
    wd00_wd13, "/Fig14_v", sbs.AssIDNo,
    "_", short.spec.lat,
    "_res", res_fac, "_mnd", r_mnd,
    ".png", sep = ""
  )
  ggsave(flnm, plot = pp,
         width = 20, height = 15, units = "cm", dpi = 300)
}
