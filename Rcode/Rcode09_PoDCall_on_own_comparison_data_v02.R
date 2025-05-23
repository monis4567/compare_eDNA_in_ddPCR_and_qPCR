#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#dev.off()
########################################################################################
# R-code for trying out the PoDCall package
# https://bioconductor.org/packages/devel/bioc/vignettes/PoDCall/inst/doc/PoDCall.html

########################################################################################
########################################################################################
#remove everything in the working environment, without a warning!!
#rm(list=ls())
########################################################################################
#install packages
#get readxl package
if(!require(readxl)){
  install.packages("readxl")
}  
library(readxl)

if(!require(readr)){
  install.packages("readr")
}  
library(readr)
#get ggplot package
if(!require(ggplot2)){
  install.packages("ggplot2")
}  
library(ggplot2)
#get magick package -as PoDCall is dependent on this package
# https://cran.r-project.org/web/packages/magick/vignettes/intro.html
if(!require(magick)){
install.packages("magick")
}  
library(magick)
#get gaston package -as PoDCall is dependent on this package
if(!require(gaston)){
  install.packages("gaston")
}  
library(gaston)

# see this website
#https://bioconductor.org/packages/devel/bioc/vignettes/PoDCall/inst/doc/PoDCall.html
# if 'PoDCall' is missing then install it
if(!require(PoDCall))
{
  ## Install from Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    library(BiocManager)
    BiocManager::install("PoDCall")
## Install PoDCall from GitHub
install.packages("devtools")
library(devtools)
# the github version requires v4.4 on 2024jan18
#devtools::install_github("HansPetterBrodal/PoDCall")
}
library(PoDCall)
#===============================================================================
# Modify the function from PoDCall -  start
#===============================================================================

checkArgumentsMultiplot <- 
  function (plateData, thresholds, channel) 
  {
    if (!is.list(plateData)) 
      stop("plateData must be a list")
    if (length(plateData) != nrow(thresholds)) 
      stop("Number of elements in\n    plateData must be equal to number of rows in thresholds")
    if (channel > 2 | channel < 1) 
      stop("invalid channel number")
    if (!is.data.frame(thresholds)) 
      stop("thresholds must be a data.frame")
    if (!("thr_target" %in% colnames(thresholds))) 
      stop("thresholds must contain\n                                                        column 'thr_target'")
    if (!("thr_ctrl" %in% colnames(thresholds))) 
      stop("thresholds must contain\n                                                        column 'thr_ctrl'")
    return(NULL)
  }

podcallMultiplot_dd <-function (plateData, thresholds, channel) 
{
  checkArgumentsMultiplot(plateData, thresholds, channel)
  if (channel == 2) {
    if (any(is.na(plateData[[1]][channel]))) {
      warning("Missing data (NA) for channel 2")
      return(NULL)
    }
    else if (any(is.na(thresholds))) {
      warning("Missing thresholds for for channel 2")
      return(NULL)
    }
  }
  plateCh <- mapply(function(x, i) {
    data.frame(wellID = i, Amplitudes = x[, c("Ch1", "Ch2")[channel]], 
               col = cut(x[, c("Ch1", "Ch2")[channel]], breaks = c(-Inf, 
                                                                   thresholds[i, c("thr_target", "thr_ctrl")[channel]], 
                                                                   Inf), labels = c("(-Inf, thr]", "[thr, Inf)")))
  }, x = plateData, i = names(plateData), SIMPLIFY = FALSE)
  plateChStacked <- rlist::list.stack(plateCh)
  thrDf <- data.frame(thresholds, wellID = names(plateCh))
  rows <- sample(nrow(plateChStacked))
  dd <- plateChStacked[rows, ]
  thrDfCh <- thrDf[, c(c("thr_target", "thr_ctrl")[channel], 
                       "wellID")]
  colnames(thrDfCh)[1] <- "thrCh"
  chCol <- c("dodgerblue3", "forestgreen")
  Amplitudes <- NULL
  wellID <- NULL
  thrCh <- NULL
  multiplot <- ggplot(data = dd, aes(x = seq_len(nrow(dd)),
                                     y = Amplitudes, group = wellID, 
                                     color = col)) + 
    geom_point(size = 1) + 
    labs(x = "Event", y = "Amplitude", color = NULL) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "none") + geom_hline(data = thrDfCh, 
                                                 aes(yintercept = thrCh), 
                                                 col = "magenta") + 
    scale_color_manual(labels = c("neg", 
                                  "pos"), values = c("gray50", chCol[channel])) + 
    facet_wrap(~wellID,  ncol = 10)
  #return(multiplot)
  lst_dd.objt<- list(dd,thrDfCh,rows,plateCh,plateChStacked)
  return(lst_dd.objt)
}



podcallMultiplot_setP <-function (plateData, 
                                  thresholds, 
                                  channel, 
                                  ncolset,
                                  nrowset,scaleset) 
{
  checkArgumentsMultiplot(plateData, thresholds, channel)
  if (channel == 2) {
    if (any(is.na(plateData[[1]][channel]))) {
      warning("Missing data (NA) for channel 2")
      return(NULL)
    }
    else if (any(is.na(thresholds))) {
      warning("Missing thresholds for for channel 2")
      return(NULL)
    }
  }
  plateCh <- mapply(function(x, i) {
    data.frame(wellID = i, Amplitudes = x[, c("Ch1", "Ch2")[channel]], 
               col = cut(x[, c("Ch1", "Ch2")[channel]], breaks = c(-Inf, 
                                                                   thresholds[i, c("thr_target", "thr_ctrl")[channel]], 
                                                                   Inf), labels = c("(-Inf, thr]", "[thr, Inf)")))
  }, x = plateData, i = names(plateData), SIMPLIFY = FALSE)
  plateChStacked <- rlist::list.stack(plateCh)
  thrDf <- data.frame(thresholds, wellID = names(plateCh))
  rows <- sample(nrow(plateChStacked))
  dd <- plateChStacked[rows, ]
  thrDfCh <- thrDf[, c(c("thr_target", "thr_ctrl")[channel], 
                       "wellID")]
  colnames(thrDfCh)[1] <- "thrCh"
  chCol <- c("dodgerblue3", "forestgreen")
  Amplitudes <- NULL
  wellID <- NULL
  thrCh <- NULL
  dd$pltrw <- substr(dd$wellID,1,1)
  dd$pltcl <- as.numeric(substr(dd$wellID,2,3))
  dd %>% dplyr::arrange(dplyr::across(pltrw,pltcl))
  multiplot <- ggplot(data = dd, aes(x = seq_len(nrow(dd)),
                                     y = Amplitudes, group = wellID, 
                                     color = col)) + 
    geom_point(size = 1) + 
    labs(x = "Event", y = "Amplitude", color = NULL) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "none") + geom_hline(data = thrDfCh, 
                                                 aes(yintercept = thrCh), 
                                                 col = "magenta") + 
    scale_color_manual(labels = c("neg", 
                                  "pos"), values = c("gray50", chCol[channel])) + 
    facet_wrap(~wellID,  ncol = ncolset, nrow=nrowset,
               scales=scaleset, dir = "h")
  return(multiplot)
  
}
#===============================================================================
# Modify the function from PoDCall -  end
#===============================================================================

# define working directory
# wd00 <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR/ddpcr_resultater"

wd00 <- getwd()
# define input file directory
wd01.1 <- "data/data_ddpcr_runs"
wd01.2 <- "amplitude_csv_data_from_QXmanager_v2"
wd00_wd01 <- paste0(wd00,"/",wd01.1,"/",wd01.2)
# define output file directory
wd02 <- "output10_figures_PoDCall"
#set working directory
#setwd(wd00)

#make complete path to input dir
wd04 <- paste(wd00,"/",wd01.1,sep="")
#make complete path to output dir
wd00_02 <- paste(wd00,"/",wd02,sep="")
# delete previous versions of the outpur directory
unlink(wd00_02, recursive=TRUE)
# # create a new version of the directory you just deleted
dir.create(wd00_02)

#_______________________________________________________________________
# section 01 -  start -  read in ddpcr QX Manager files
#_______________________________________________________________________
# define prefix for  w input directories
prefix_ind_01 = "dd"
#list all directories in input directory 
ls.d01 <- list.dirs(path = wd00_wd01, full.names = TRUE, recursive = TRUE)
# get a list of directories that has the postfix for input files
ind_01 <- ls.d01[grep(prefix_ind_01,ls.d01)]

# species abbreviation to isolate for 
spc.Abbr_iso <- "Mnelei"
spc.Abbr_iso <- "Erisin"

# make a vector of the species abbreviations to iterate over
set_of_spc.Abbr <- c("Mnelei","Erisin","Bonham","Calsap","Promin","Psever")
#set_of_spc.Abbr <- c("Psever", "Psefar")
#set_of_spc.Abbr <- c("Psefar")
# get the total number of elements in the vector containing species 
# abbreviations names
nspcAbr <- length(set_of_spc.Abbr)
#nspcAbr <- 1
# make an empty list that plots and data frame can be collected in
lst_plts <- list()
lst_PDCe_dds  <- list()
# iterate over elements in the vector of numbers for species abbreviations
for (spcAbbr.no in seq(1,nspcAbr,1))
{
  # get the species abbreviation that matches the species number
  spc.Abbr_iso <- set_of_spc.Abbr[spcAbbr.no]
# then isolate for only this species abbreviation
ind_01 <- ls.d01[grep(spc.Abbr_iso,ls.d01)]
#ind_01 <- ls.d01[grep("Psefar",ls.d01)]
ind_01 <- ind_01[1]
# define the directory with the input files
# note that 'podcallDdpcr' requires a directory that is full of
# csv files, where each individual csv files represents one single well
# in the plate. Such amplitude csv files can only be exported by
# QX manager v2 or higher. I only had access to the QX manager v1.2
# to begin with, and was therefore without the correct csv amplitude files
# on
# ## Run PoDCall
# thresholdTable <- PoDCall::podcallDdpcr(dataDirectory=ind_01, 
#                                          software="QX Manager",
#                                          ch2=T)


#4.2 importSampleSheet()
#Reads a .csv-file outputted from QuantaSoft or QX Manager
#to get information about the samples: Sample name/id, Assay for
#target and control.


## Select wells to get information for
# define directory with csv files with sample sheets from
# the QX Manager v1.2 software. Note that versions >2 of the QX Manager
# also can be used. I am just unfortunate and only have access to v1.2
wd03 <- "csv_sample_sheets_from_QXmanager"
# paste work directory and path for directory for the csv files with 
# sample sheets
wd00_03 <- paste0(wd00_wd01,"/",wd03)
#make all combinations of well names
dfwells <- expand.grid(seq(1:7),LETTERS[1:8])
#pad with zeros to 2 characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
dfwells$Var1 <- stringr::str_pad(dfwells$Var1, 2, pad = "0")
w.id <- paste0(dfwells$Var2,dfwells$Var1)

# list all files
ls.ssht <- list.files(wd00_03)
ssht <-  ls.ssht[grep(paste0(spc.Abbr_iso,"_"),ls.ssht)]
p.ssht <- paste0(wd00_03,"/",ssht)
## Read in sample sheet information for selected wells
ssht01 <- PoDCall::importSampleSheet(sampleSheet=p.ssht, 
                                     well_id=w.id,
                                     software="QuantaSoft")
                                 #software="QX Manager")

# list all files with xls setup for ddPCR
ls.ssht <- list.files(wd04)
xlssht <-  ls.ssht[grep("\\.xls",ls.ssht)]
xlssht <-  ls.ssht[grep("plate_setup_ddpcr",ls.ssht)]
xlssht <-  paste0(wd00,"/",xlssht)

# get the specific sample sheet
xlssht.s <- xlssht[grep(spc.Abbr_iso,xlssht)]

library(readxl)
spsht.stup <- readxl::read_excel(xlssht.s)
# get the row number where the setup of well starts
rwn.wtbs <- which("WellNumber"==spsht.stup[,1])
# get a seq from where 'Well' appears and the next 96 rows
setuprws <- seq(rwn.wtbs,(rwn.wtbs+96))
# limit the setup file to only comprise these rows with the setup
ddstup <- spsht.stup[setuprws,]
# make the setup data frame a data frame instead of a tibble 
df_dds <- as.data.frame(ddstup)
# use the 1st row as column names 
colnames(df_dds) <- df_dds[1,]
# get the data frame without the 1st row
df_dds <- df_dds[-1,]
# make a column that has 'WellNo' as column name, in the qpcr setup data frame
df_dds$WellNo <- df_dds$WellNumber
# get the sampl_id from the excel spreadsheet data frame that has the ddPCR
# setup table
ssht01$sample_id <- df_dds$WellName[match(ssht01$well_id,df_dds$WellNo)]
ssht01$target_assay <- spc.Abbr_iso
ssht01$ctrl_assay <- "None"
# subset the data frame to only comprise the standard dilution series
ssht01.stds <- ssht01[grepl("std",ssht01$sample_id),]
ssht01.stds$stdllvl <- gsub(".*std_(.*)","\\1",ssht01.stds$sample_id)

#_______________________________________________________________________________


# define the directory that has the sub directories with amplitude data
wd_ampl <- "amplitude_csv_data_from_QXmanager_v2"
# paste together with the working directory
wd00_ampl <- paste0(wd00_wd01,"/",wd_ampl)
lstdrs.wd_ampl <- list.dirs(path =wd00_ampl, full.names = TRUE, recursive = TRUE)
# get the directory that has the amplitude data for the species
dir.A_sp  <- lstdrs.wd_ampl[grepl(spc.Abbr_iso,lstdrs.wd_ampl)]

# 4.3 - here : https://bioconductor.org/packages/devel/bioc/vignettes/PoDCall/inst/doc/PoDCall.html
## Path to example data files included in PoDCall

#path <- system.file("extdata", "Amplitudes/", package="PoDCall")
## Read in data files
# Notice that the amplitude data file begins with
# Target Value of 0 = negative
# Target Value of 1 = positive
# Target Value of u = unclassified (Advanced Classification Mode)
# 
# Ch1Amplitude,Ch2Amplitude,Mnelei
# 4219.069,314.6323,1
# Which mean you will need to skip the first 4 lines
ampData <- importAmplitudeData(dataDirectory=dir.A_sp, skipLines=4)
# for some reason some of the individual data frames in th elong list of data
# frames ende uo with columns that has characters insterad of numeric values
# this makes the next step with 'podcallThresholds' going into error
# the solution was to ensure all columns in every data frame in the list
# of data frames all have numeric values instead of characters
ampData <- lapply(ampData, dplyr::mutate_if, is.character, as.numeric)

## Now that the amp data is ensured to have only numeric values, the next 
## step is to 
## Calculate thresholds, metrics, concentrations
thresholdTable <- PoDCall::podcallThresholds(plateData=ampData, 
                                     #nchannels=c(1,2)[2],
                                    B=200, 
                                    Q=9,
                                    refWell = 1)

#print(thresholdTable)
thrTable <- thresholdTable

#4.4 -podcallChannelPlot() 
## Read in data and threshold table
## path <- system.file("extdata", "Amplitudes/", package="PoDCall")
# ampData <- importAmplitudeData(dataDirectory=dir.A_sp, skipLines=4)


# #data("thrTable")
# thresholdTable <- thrTable
## Select channel and well to plot
ch <- 1 # target channel
well_id <- names(ampData)[1] # First well in list
## Plot title
plotTitle <- paste0(well_id, ", Ch", ch)
# ## Create plot
# podcallChannelPlot(channelData=ampData[[well_id]][,ch],
#                    thr=thresholdTable[well_id, "thr_target"],
#                    channel=ch,
#                    plotId=plotTitle)


#4.5 podcallScatterplot()
## Read in data and threshold table
#path <- system.file("extdata", "Amplitudes/", package="PoDCall")
# ampData <- importAmplitudeData(dataDirectory=dir.A_sp, skipLines=4)
# thresholdTable <- thrTable
## Select channel and well to plot
ch <- 1 # target channel
well_id <- names(ampData)[1] # First well in list
## Plot title
plotTitle <- paste0(well_id, ", Ch", ch)
## Create plot
# podcallScatterplot(channelData=ampData[[well_id]][,ch],
#                    thr=thresholdTable[well_id, "thr_target"],
#                    channel=ch,
#                    plotId=plotTitle)




# 4.6  podcallHistogram()
## Read in data and threshold table
#path <- system.file("extdata", "Amplitudes/", package="PoDCall")
# ampData <- importAmplitudeData(dataDirectory=dir.A_sp, skipLines=4)
# thresholdTable <- thrTable
## Select channel and well to plot
ch <- 1 # target channel
well_id <- names(ampData)[1] # First well in list
## Plot title
plotTitle <- paste0(well_id, ", Ch", ch)
## Create plot
# podcallHistogram(channelData=ampData[[well_id]][,ch],
#                  thr=thresholdTable[well_id, "thr_target"],
#                  channel=ch,
#                  plotId=plotTitle)


# unique(is.infinite(unlist(ampData)))
#4.7  podcallMultiplot()
## Read in data and threshold table
#path <- system.file("extdata", "Amplitudes/", package="PoDCall")
# ampData <- importAmplitudeData(dataDirectory=dir.A_sp, skipLines=4)
# thresholdTable <- thrTable
## Channel to plot
ch <- 1

# ## Create comparison plot
# podcallMultiplot(plateData=ampData,
#                  thresholds=thresholdTable[names(ampData),], 
#                  channel=ch)
#_______________________________________________________________________________
## Create comparison plot elements but using the modified function, that
# makes a list of the components used for the plotting function
PDCe <- podcallMultiplot_dd(plateData=ampData,
                 thresholds=thresholdTable[names(ampData),], 
                 channel=ch)

chCol <- c("dodgerblue3", "forestgreen")
channel <- ch
# use the modified version of the function
# to get the individual elements
PDCe_dd<- as.data.frame(PDCe[[1]])
PDCe_dd %>% dplyr::arrange(wellID)
PDCe_thrDfCh <- as.data.frame(PDCe[[2]])
PDCe_dpcnt <- as.data.frame(PDCe[[3]])
PDCe_plateCh <- (PDCe[4])
PDCe_plateChStacked <- as.data.frame(PDCe[[5]])


# subset among these individual elements to only get 
# samples that were analysed in the '01' column of the plate
PDCe_dd.01 <- PDCe_dd[grepl("01",PDCe_dd$wellID),]
thrDfCh.01 <- PDCe_thrDfCh[grepl("01",PDCe_thrDfCh$wellID),]
# mathc to get the standard dilution level
PDCe_dd.01$stdllvl <- ssht01.stds$stdllvl[match(PDCe_dd.01$wellID, ssht01.stds$well_id)]
# assume that all NAs are NegControls in this column 
PDCe_dd.01$stdllvl[is.na(PDCe_dd.01$stdllvl)] <- 'NegC'


# Now make a plot only based on this subsetted version of the data frame
# now you can control the gplot2 elements yourself
plt02 <- ggplot(data = PDCe_dd.01) + 
  theme_bw() +
  geom_point(aes(x = seq_len(nrow(PDCe_dd.01)),
                 y = Amplitudes, group = wellID, 
                 color = col),size = 1) + 
  labs(x = "Event", y = "Amplitude", color = NULL) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none") + 
  geom_hline(data = thrDfCh.01, 
              aes(yintercept = thrCh), 
              col = "magenta") + 
  scale_color_manual(labels = c("neg", 
                                "pos"), 
                     values = c("gray50", 
                                chCol[channel])) +
  # Add a label that can hold the dilution step information
  # https://stackoverflow.com/questions/36089962/ggplot2-change-strip-text-position-in-facet-grid-plot
  geom_text(aes(label = stdllvl),
            col="black",
            x = Inf, y = Inf, 
            hjust = 1.5, vjust = 1.5, check_overlap = TRUE) +
  # reduce the spacing between the panels in facet wrap
  theme(panel.spacing = unit(0.1, "lines")) +
  facet_wrap(~wellID,  ncol = 8) +
  ggtitle(spc.Abbr_iso)
  
#plt02
# place the plot in the list that collects the plots
lst_plts[[spcAbbr.no]] <- plt02
# get the 'thrCh' value per well ID and add it to the data frame
PDCe_dd.01$thrCh <- thrDfCh.01$thrCh[match(PDCe_dd.01$wellID,thrDfCh.01$wellID)]
# also add the species abbreviation to the data frame
PDCe_dd.01$spc.Abbr  <- spc.Abbr_iso
# also add the chCol[channel] to the data frame
PDCe_dd.01$chCol  <- chCol[channel]
# finally collect the prepared data frames in a list
lst_PDCe_dds[[spcAbbr.no]] <- PDCe_dd.01
# end the iteration over the species abbreviations
}


#bind the rows in each list in to one data frame
df_PDCes <- data.table::rbindlist(lst_PDCe_dds, fill=T)
df_PDCes <- as.data.frame(df_PDCes)

#_______________________________________________________________________________
# best solution for multiple plots no 1 -start
#_______________________________________________________________________________
# see : https://stackoverflow.com/questions/62652308/combine-multiple-facet-strips-across-columns-in-ggplot2-facet-wrap

library(tidyverse)
library(gtable)
library(grid)

idx = 1:16

p1 = expand_grid(id=idx, id2=c("A", "B"), x=1:10) %>%
  mutate(y=rnorm(n=n())) %>%
  ggplot(aes(x=x,y=y)) +
  geom_jitter() +
  facet_wrap(~id + id2, nrow = 4, ncol=8)

g <- ggplot_gtable(ggplot_build(p1))

stript <- grep("strip", g$layout$name)

grid_cols <- sort(unique(g$layout[stript,]$l))
t_vals <- rep(sort(unique(g$layout[stript,]$t)), each = length(grid_cols)/2)
l_vals <- rep(grid_cols[seq_along(grid_cols) %% 2 == 1], length = length(t_vals))
r_vals <- rep(grid_cols[seq_along(grid_cols) %% 2 == 0], length = length(t_vals))
labs   <- levels(as.factor(p1$data$id))

for(i in seq_along(labs))
{
  filler <- rectGrob(y = 0.7, height = 0.6, gp = gpar(fill = "gray80", col = NA))
  tg    <- textGrob(label = labs[i], y = 0.75, gp = gpar(cex = 0.8))
  g     <- gtable_add_grob(g, filler, t = t_vals[i], l = l_vals[i], r = r_vals[i], 
                           name = paste0("filler", i))
  g     <- gtable_add_grob(g, tg, t = t_vals[i], l = l_vals[i], r = r_vals[i], 
                           name = paste0("textlab", i))
}

grid.newpage()
grid.draw(g)
#_______________________________________________________________________________
# best solution for multiple plots no 1 -end
#_______________________________________________________________________________

library(tidyverse)
library(gtable)
library(grid)
library(patchwork)
#define working directory
#wd00.1 ="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6"
wd00.1 ="/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR/MONIS6_2021_data"

#define input directory
wd01.1 <- "output02_merged_txtfiles_from_mxpro_for_MONIS6"
wd02.1 <- paste(wd00.1,"/",wd01.1,sep="")
#define inputfile
inpf01 <- "outfile02_merged_mxpro_csvfls_MONIS6.csv"
#paste path and input file together
if01wp <- paste(wd02.1,"/",inpf01,sep="")
# read in an xlsx file with all primer assays listed 
# for each species anmd the abbreviations
inf04 <- "list_of_specific_assays_MONIS6.xlsx"
wd00.1_inf04 <- paste0(wd00.1,"/",inf04 )
library(xlsx)
df_dtc_asss <- xlsx::read.xlsx(wd00.1_inf04,1)

# length Latin Genus and Species names
df_PDCes$GnNm <-  df_dtc_asss$Genus[match(df_PDCes$spc.Abbr,df_dtc_asss$AbbrvNm)]
df_PDCes$SpNm <-  df_dtc_asss$Species[match(df_PDCes$spc.Abbr,df_dtc_asss$AbbrvNm)]
df_PDCes$GnNmLt <- substr(df_PDCes$GnNm, 1, 1)
df_PDCes$speciesNm <- paste0(df_PDCes$GnNmLt,". ",df_PDCes$SpNm)
#unique(df_PDCes$speciesNm)

# use the multiple plots setup on the actual data
#  PDCe_dd.01$chCol 
p2 <- ggplot(data = df_PDCes) + 
  # # justify to the left margin
  theme_minimal() +
  # 
  #   #https://stackoverflow.com/questions/62009919/facet-title-alignment-using-facet-wrap-in-ggplot2
  #   theme(strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
  #         panel.background = element_rect(fill = "grey99",
  #                                         color = "white")) +
  
  geom_point(aes(x = seq_len(nrow(df_PDCes)),
                 y = Amplitudes, group = wellID, 
                 color = col),size = 1) + 
  labs(x = "Event", y = "Amplitude", color = NULL) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  geom_hline(#data = thrDfCh.01, 
    aes(yintercept = thrCh), 
    col = "magenta") + 
  scale_color_manual(labels = c("neg", 
                                "pos"), 
                     values = c("gray50", 
                                chCol[channel])) +
  # Add a label that can hold the dilution step information
  # https://stackoverflow.com/questions/36089962/ggplot2-change-strip-text-position-in-facet-grid-plot
  geom_text(aes(label = stdllvl),
            col="black",
            x = Inf, y = Inf, 
            hjust = 1.5, vjust = 1.5, check_overlap = TRUE) +
  # reduce the spacing between the panels in facet wrap
  theme(panel.spacing = unit(0.1, "lines")) +
  #facet_wrap(~ spc.Abbr + wellID,  nrow = 4, ncol=8) +
  facet_wrap(~ spc.Abbr + wellID,  nrow = nspcAbr, ncol=8) +
  #ggtitle("nice plot")
  ggtitle(" ")


g <- ggplot_gtable(ggplot_build(p2))

stript <- grep("strip", g$layout$name)

grid_cols <- sort(unique(g$layout[stript,]$l))
t_vals <- rep(sort(unique(g$layout[stript,]$t)), 
              each = length(grid_cols)/8)
l_vals <- rep(grid_cols[seq_along(grid_cols) %% 8 == 1], 
              length = length(t_vals))
r_vals <- rep(grid_cols[seq_along(grid_cols) %% 8 == 0], 
              length = length(t_vals))
labs   <- levels(as.factor(p2$data$speciesNm))

for(i in seq_along(labs))
{
  filler <- rectGrob(y = 0.7, height = 0.6, just = "centre",
                     gp = gpar(fill = "grey76", col = NA))
  tg    <- textGrob(label = labs[i], y = 0.75, just = "left",
                    gp = gpar(cex = 0.8))
  g     <- gtable_add_grob(g, filler, 
                           t = t_vals[i], l = l_vals[i], r = r_vals[i], 
                           name = paste0("filler", i))
  g     <- gtable_add_grob(g, tg, 
                           t = t_vals[i], l = l_vals[i], r = r_vals[i], 
                           name = paste0("textlab", i))
}

# grid.newpage()
# p3 <- grid.draw(g)
bSaveFigures<-T
if(bSaveFigures==T){
  ggsave(plot = g, 
         filename = paste0(wd00_02,"/Fig01_v01_PoDCall_plot_std_dil_ser.png"),
         #width=210*0.64,height=297,
         #width=210*0.8,height=297,
         #width=297,height=210,
         width=210,height=297,
         #width=1.4*297,height=210,
         units="mm",dpi=300)
}

#____________

# Check out this webpage on tips for gathering plots from ggplot
# https://r-charts.com/ggplot2/combining-plots/#gridExtra
# https://wilkelab.org/cowplot/articles/plot_grid.html

# install.packages("ggplot2")
# install.packages("cowplot")
library(ggplot2)
library(cowplot)
# 
# plot_grid(lst_plts[[1]], 
#           lst_plts[2],
#           lst_plts[3],
#           labels = c('A', 'B','C'), # Or "AUTO" or "auto"
#           label_fontfamily = "serif",
#           label_fontface = "bold",
#           label_colour = "dodgerblue2")


