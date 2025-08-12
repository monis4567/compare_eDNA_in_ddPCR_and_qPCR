#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
# “comparison of qPCR and ddPCR”
# Authors: Steen Wilhelm Knudsen.
# Change the working directory to a path on your own computer , and run
# the individual parts below to reproduce the diagrams presented in the paper
#
# All input data required needs to be available as csv-files in the same directory 
# as this R-code use for working directory.
#
# Occassionally the code will have difficulties producing the correct diagrams,
# if the packages and libraries are not installed.
# Make sure the packages are installed, and libraries are loaded, if the R-code
# fails in producing the diagrams.
#
#________________IMPORTANT!!_________________________________________________#
# (1)
#You have to change the path to the working directory before running this code
#
# (2)
# The 4 data input files required:
#
# must be located in the same working directory - as specified in the code below
         
#____________________________________________________________________________#
#____________________________________________________________________________#
# R-code provided for the project:
#remove everything in the working environment, without a warning!!
#rm(list=ls())
#get rinat package
# if(!require(rinat)){
#   remotes::install_github("ropensci/rinat")
#   install.packages("rinat")
# }  
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(rinat)
#define working directory
# wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2022/ddPCR_qPCR_MST"
# wd00  <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR"
wd00   <- getwd()
#define input file  directory
wd01 <- "data/data_ddpcr_runs"
# paste together working directory and input file  directory
wd00_01 <- paste0(wd00,"/",wd01)

# make an output directory
wdout <- "output12_map_iNat_records_and_ddpcr"
#make complete path to output dir
wd00_wdout <- paste(wd00,"/",wdout,sep="")
#Delete any previous versions of the output directory
unlink(wd00_wdout, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wdout)



#_______________________________________________________________________________
# 01 start - make trycatch function, in case the genus name is not on GBIF 
#_______________________________________________________________________________
# https://www.statology.org/r-trycatch/
try.c.iNat <- function(tx, boundslim, endval){
  tryCatch(
    {
      #get iNaturalist records
      g <- rinat::get_inat_obs(
        taxon_name = tx,
        quality = "research",
        geo=T, #only include geo referenced results
        bounds = boundslim,
        maxresults = endval)
      return(g)
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
}
#_______________________________________________________________________________
# 01 end - make trycatch function, in case the genus name is not on GBIF 
#_______________________________________________________________________________


# https://github.com/cran/rinat
#_______________________________________________________________________________
# nelng cuts on the eastern boundary
# nelat cut on the northern border
#                             |
#                       nelng |
#             N         nelat___    |
#             |                     | y-axis is lat
#             |                     |
#             |                     |
#   W____________________E          |
#             |                     |
#             |                     |
#             |                     |
# ___  swlat  S                     |
#     |swlng
#     |
#_____________________________
#             x-axis is lon
# swlng cuts on the western boundary
# swlat cut on the southern border
#try defining your own bounding box
set_nelat= 58
set_nelng= 15.4
set_swlat= 54.4
set_swlng= 8
#try defining your own bounding box
set_nelat= 59.5
set_nelng= 17
set_swlat= 53
set_swlng= 7
# needs to be in the format 
#  'min_lat','min_lon','max_lat', 'max_lon'
# whcih equals
#  'min_y','min_x','max_y', 'max_x'
boundslim <- c(set_swlat, set_swlng, set_nelat, set_nelng)
## Search using the boundslimits defined above
df_iNat01 <- get_inat_obs(taxon_name = "Mya arenaria",
                          quality = "research",
                          bounds = boundslim,
                          maxresults = 500)
# make the ggplot
plt_amp01 <- ggplot(data = df_iNat01, aes(x = longitude,
                                          y = latitude,
                                          colour = scientific_name)) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "grey95",
               color = "gray40",
               linewidth = 0.1) +
  geom_point(size = 0.7, alpha = 0.5) +
  coord_fixed(xlim = range(df_iNat01$longitude, na.rm = TRUE),
              ylim = range(df_iNat01$latitude, na.rm = TRUE)) +
  theme_bw()
# see the plot
plt_amp01



#_______________________________________________________________________________

#make a list of organisms
lst_orgnsm <- c(
  "Mya arenaria",
  "Callinectes sapidus")
#_______________________________________________________________________________
# Get the name for data frame as
# written out as a csv file from the previous R code
wd11 <- "output11_compare_all_ddpcr_and_qpcr"
# make a file name to read in
filNm.e09 <- paste0(wd00,"/",wd11,"/table_w_ddpcr_and_qpcr_data_and_lat_lon_v01.csv")
# read in the saved  data frame
dfg_dq <- read.csv(file=filNm.e09)
#paste together genmus name and species name
dfg_dq$Lat_Species <- paste0(dfg_dq$GnNm," ",dfg_dq$SpNm)
# vmake a alist of the species names
latSpcNms  <- unique(dfg_dq$Lat_Species)
# limit to species names that have a 'sapce' included
latSpcNms <- latSpcNms[(grepl(" ",latSpcNms))]
#_______________________________________________________________________________
# copy the list into a new object
lst_orgnsm <- latSpcNms

lst_spcs <- lst_orgnsm
#_______________________________________________________________________________

# 01 start - make trycatch function, in case the genus name is not on GBIF 
#_______________________________________________________________________________
# https://www.statology.org/r-trycatch/
try.c.iNat <- function(tx, boundslim, endval){
  tryCatch(
    {
      #get iNaturalist records
      g <- rinat::get_inat_obs(
        taxon_name = tx,
        quality = "research",
        geo=T, #only include geo referenced results
        bounds = boundslim,
        maxresults = endval)
      return(g)
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
}
#_________________________
# 01 end - make trycatch function, in case the genus name is not on GBIF 
#_________________________
kee <- c("scientific_name",
         "datetime",
         "place_guess",
         "latitude",
         "longitude",
         "tag_list",
         "common_name",
         "url",
         "image_url",
         "species_guess",
         "iconic_taxon_name",
         "taxon_id",
         "num_identification_agreements",
         "num_identification_disagreements",
         "observed_on_string",
         "observed_on",
         "time_observed_at",
         "time_zone")

#make an empty list to use for collecting data frame
lst_tx_gobs <- list()
#start a growing number
i <- 1
lst_spcs_a <- lst_spcs[3]
lst_spcs_a <- lst_spcs
#iterate over taxon names in list 
for (tx in lst_spcs_a)
{
  #  print(tx)}
  print(tx)
  # substitute the underscore
  tx <- gsub("_"," ",tx)
  ## Search for the  species using the boundslimits defined above
  g <- try.c.iNat(tx, boundslim, 5000)
  # limit to only specific columns, otherwise it ends up being too much data
  g <- g[kee]
  # check if there are no data, and in that case, add NAs for the 'kee' columns
  if(!is.null(colnames(g)))
  {g <- g} else {
    df_tmp <- as.data.frame(t(as.matrix(kee)))
    df_tmp <- rbind(df_tmp,rep(NA,length(kee)))
    colnames(df_tmp) <- df_tmp[1,]
    df_tmp <- df_tmp[-1,]
    g <- df_tmp
  }
  # add the taxon name that was used for making the search
  g$txNmsrch <- tx
  # make the entire data frame characters
  g[] <- lapply(g, as.character)
  # append the data frame to the list of data frames
  # store it as the i'th element 
  lst_tx_gobs[[i]] <- g
  # increase the count of i by one
  i <- i+1
  # end iteration over amphibian species in the list
}
#bind the rows in each list in to one data frame
df_g03 <- data.table::rbindlist(lst_tx_gobs, fill=T)
df_g03 <- as.data.frame(df_g03)
# if there is no latitude, then omit the row
df_g03 <- df_g03[!is.na(df_g03$lat),]
# copy the column with the scientific name
df_g03$Lat_Species <- df_g03$scientific_name

#_______________________________________________________________________________
# define columns to keep
ck.03 <- c( "latitude",
            "longitude",
            "observed_on",
            "Lat_Species")
# keep only columns listed in vector
df_g04 <- df_g03[ck.03]
# re order the columns by alphabetical order
df_g04 <- df_g04[order(colnames(df_g04))]
# add a column that defines that this is an iNaturalist observation
df_g04$source <- "iNat"
df_g04$eval_detect <- "iNat"
# make an evaluation column , that by default sets 
# the value to 'no eDNA' detected
dfg_dq$eval_detect <- "no eDNA"
# evaluate if eDNA was detected - i.e. the 'eDN.lv1.l10' is not zero 
dfg_dq$eval_detect[(dfg_dq$eDN.lv1.l10!=0)] <- "eDNA detected" 

# define columns to keep
c_dq.k <- c( "lon",
             "lat",
             "Dato_inds",
             "mch",
             "Lat_Species",
             "eval_detect")
# keep only columns listed in vector
dfdq_02<- dfg_dq[c_dq.k]


# rename column names
colnames(dfdq_02) <-  c(
            "longitude",
            "latitude",
            "observed_on",
            "source",
            "Lat_Species",
            "eval_detect")
# re order the columns by alphabetical order
dfdq_02 <- dfdq_02[order(colnames(dfdq_02))]
# get the column names for the data frames
clNmd04 <-colnames(df_g04)
clNmd02 <- colnames(dfdq_02)
# compare the columns names
clNmd04 %in% clNmd02
# # make a column that denotes if eDNA was detected
# dfdq_02$eval_detect <- "eDNA_present"
# bind rows onto the data frame to have the two data frames combined
df_A06 <- rbind(dfdq_02,df_g04)
#View(df_A06)
# check the evaluation categories in the data frame
unique(df_A06$eval_detect)
#_______________________________________________________________________________
library(dplyr)
# get the number of latin species
LSpc <- unique(dfg_dq$Lat_Species)

# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# install.packages("rnaturalearthhires")
# remotes::install_github("ropensci/rnaturalearthhires")
# install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")

library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

# split the string
obsdtspl <- strsplit(as.character(df_A06$observed_on), "-")
# and use the splitted string to get each element
# and make numeric
obd.y <- as.numeric(sapply(obsdtspl, "[[", 1))
obd.m <- as.numeric(sapply(obsdtspl, "[[", 2))
obd.d <- as.numeric(sapply(obsdtspl, "[[", 3))
# combine to a data frame
df_obs <- as.data.frame(cbind(obd.y,
                              obd.m,
                              obd.d))
# add the date for year, month and day of observation
# by column binding the data frames
df_A06 <- cbind(df_obs,df_A06)
# make an empty column for season evaluations
df_A06$ssnno <- NA
# evalutate on the season column to add a number version of the season
df_A06$ssnno[(df_A06$obd.m<=6)] <- "1st"
df_A06$ssnno[(df_A06$obd.m>6)] <- "2nd"
# substitute the season number with a string
df_A06$ssnno2 <- gsub("1st","Jan-Jun",df_A06$ssnno)
df_A06$ssnno2 <- gsub("2nd","Jul-Nov",df_A06$ssnno2)
# copy the column
df_A06$yer_ssn2 <- df_A06$ssnno2
# make sure the latitude and longitude are numeric
df_A06$latitude <- as.numeric(df_A06$latitude)
df_A06$longitude <- as.numeric(df_A06$longitude)
# 
unique(dfg_dq$mch.ssn)

dfg_dq %>% dplyr::select()
# make a sequence for this range
nfLSpc<- seq(1,length(LSpc),1)
#nfLSpc <- 13
# iterate over species
for (i in nfLSpc)
{print(paste0(i, " making plot for ",LSpc[i]))
  #}
  # i <- 13
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  ins <- stringr::str_pad(i, 2, pad = "0")
  # get the species name with an underscore
  sbs_spcNm_wu <- gsub(" ","_",LSpc[i])
  #}
  # make a function to use for making facet wrap headers
  # https://ggplot2.tidyverse.org/reference/as_labeller.html
  # https://stackoverflow.com/questions/63857833/changing-the-facet-wrap-labels-using-labeller-in-ggplot2
  #appender <- function(nm1) paste0(df_tx01$class[match(nm1,df_tx01$family )],": ", nm1)
  
  df_A06.1 <- df_A06[grepl(LSpc[i],df_A06$Lat_Species),]
  # copy the data frame
  df_A06.2 <-  df_A06.1
    # add a level for jittering points
  jitlvl <- 0.03
  # subset data frame
  df_A07 <- df_A06.2[(df_A06.2$eval=="iNat"),]
  
  # make a data frame 
  df_A07.ne <- df_A06.1[(df_A06.2$eval=="no eDNA"),]
  df_A07.ed <- df_A06.1[(df_A06.2$eval=="eDNA detected"),]
  
  #make plot
  p05 <- ggplot(data = world) +
    geom_sf(color = "black", fill = "azure3", lwd=0.1) +
    # also see: https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
    theme_void() +
    #https://ggplot2.tidyverse.org/reference/position_jitter.html
    #https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2
    # use 'geom_jitter' instead of 'geom_point' 
    geom_jitter(data = df_A07 ,
                aes(x = longitude, y = latitude),
                color=alpha(c("black"),c(0.6)),
                size=1.6,
                fill=alpha(c("yellow"),c(0.6)),
                shape=22,
                width = jitlvl, #0.07, jitter width 
                height = jitlvl) + #, #0.07, # jitter height
    #Arrange in facets
    ggplot2::facet_wrap( ~ yer_ssn2,
                         drop=FALSE,
                         dir="h",
                         ncol = 2,
                         labeller = label_bquote(cols =
                                                   .(as.character(yer_ssn2))
                         ) ) +
    # alter the them strip above
    theme(strip.text = element_text(#face = "bold",
      color = "black",
      hjust = 0,
      size = 8),
      strip.background = element_rect(fill = c("white"),
                                      #linetype = "solid",
                                      color = "white",
                                      linewidth = 1)) +
    # ggplot2::facet_wrap( ~ Lat_Species,
    #                      drop=FALSE,
    #                      ncol = 4,
    #                      labeller = as_labeller(appender))+
    #define limits of the plot 
    ggplot2::coord_sf(xlim = c(6.6, 17.2),
                      ylim = c(54.2, 58.4), 
                      expand = FALSE) +
    # # Add points for sampling locations -  this has to be from a separate data frame
    # The 'shape=3' makes crosses for sampled locations without eDNA
    geom_point(data=df_A07.ne,aes(x=longitude,y=latitude),
               shape=3,colour=alpha("#000000",1),size=1.6) +
    # The 'shape=21' makes circular points for sampled locations with eDNA present
    geom_point(data=df_A07.ed,aes(x=longitude,y=latitude),
               shape=21,colour=alpha("#000000",1),fill=alpha("firebrick3",0.6),
               size=2.4) +
    labs(title = LSpc[i]) +
    # make the  labels turn 90 degrees on the x-axis
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # or remove them all completely
    #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()) +
    # legend on bottom
    theme(legend.position="bottom")  + 
    # make the title in italic
    theme(plot.title = element_text(face = "italic", size =10))  
  
  #change axis labels
  p05t <- p05 + xlab("longitude") + ylab("latitude")
  #p05t <- p05 + xlab("længdegrad") + ylab("breddegrad")
  #p05t <- p05 + xlab(" ") + ylab(" ")
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  # p05t <- p05t + labs(color='')
  # p05t <- p05t + labs(fill='')
  # p05t <- p05t + labs(shape='')
  # p05t <- p05t + labs(size='')
  #get the number of species
  ncat <- length(unique(df_A07$eval))
  # https://github.com/tidyverse/ggplot2/issues/3492
  #repeat 'black' a specified number of times
  filltxc = rep("black", ncat)
  #adjust tick marks on axis
  p05t <- p05t + scale_y_continuous(breaks=seq(54.2, 58.4,2))
  p05t <- p05t + scale_x_continuous(breaks=seq(6.6, 17.2,4))
  
  # see the plot
  p05t
  
  #
  bSaveFigures<-T
  if(bSaveFigures==T){
    ggsave(plot = p05t, 
           filename = paste0(wd00_wdout,"/Fig11_v",ins,"_map_of_",sbs_spcNm_wu,"_detected_2017_to_2023.png"),
           width=210*0.60,height=297*0.17,
           #width=210*0.8,height=297,
           #width=297,height=210,
           #width=297,height=210,
           #width=1.4*297,height=210,
           units="mm",dpi=300)
  }
  # end iteration over latin species names
}

# write the data frame to a csv file
write.csv(df_A06,paste0(wd00_wdout,
                        "/table_w_ddpcr_results_v02.csv"))


