#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

########################################################################################
# R-code for trying out the PoDCall package
# https://bioconductor.org/packages/devel/bioc/vignettes/PoDCall/inst/doc/PoDCall.html

########################################################################################
########################################################################################
#remove everything in the working environment, without a warning!!
rm(list=ls())
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

# define working directory
wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2022/ddPCR_qPCR_MST/ddpcr_resultater"
wd00 <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR"
wd00 <-  paste0(wd00,"/ddpcr_resultater")
# define input file directory
wd01 <- "amplitude_csv_data_from_QXmanager_v2"
# define output file directory
wd02 <- "output02_figures_PoDCall"
#set working directory
setwd(wd00)

#make complete path to input dir
wd00_01 <- paste(wd00,"/",wd01,sep="")
#make complete path to output dir
wd00_02 <- paste(wd00,"/",wd02,sep="")
# delete previous versions of the outpur directory
unlink(wd00_02, recursive=TRUE)
# # create a new version of the directory you just deleted
dir.create(wd00_02)


# 3
# This file is the common traditional csv outputfile from the QX
# Manager software - bothe QX Manager v1.2. and v2.1 are able
# to produce this csv file
## Path to example files included in PoDCall
path <- system.file("extdata", "Sample_names.csv", package="PoDCall")
Amplt.path <- "/home/hal9000/R/x86_64-pc-linux-gnu-library/4.3/PoDCall/extdata/Amplitudes"
## Run PoDCall
thresholdTable <- PoDCall::podcallDdpcr(dataDirectory=Amplt.path, 
                               software="QuantaSoft")

# 4.1 - here : https://bioconductor.org/packages/devel/bioc/vignettes/PoDCall/inst/doc/PoDCall.html
## Path to example data files included in PoDCall
path <- system.file("extdata", "Amplitudes/", package="PoDCall")
## Read in data files
dataList <- PoDCall::importAmplitudeData(dataDirectory=path, skipLines=0)
str(dataList)

# 4.2 - here : https://bioconductor.org/packages/devel/bioc/vignettes/PoDCall/inst/doc/PoDCall.html
## Path to example files included in PoDCall
path <- system.file("extdata", "Sample_names.csv", package="PoDCall")
## Select wells to get information for
well_id <- c("A04", "B04", "D04")
## Read in sample sheet information for selected wells
sampleSheet <- PoDCall::importSampleSheet(sampleSheet=path, well_id=well_id,
                                 software="QuantaSoft")
print(sampleSheet)

# 4.3 - here : https://bioconductor.org/packages/devel/bioc/vignettes/PoDCall/inst/doc/PoDCall.html
## Path to example data files included in PoDCall
path <- system.file("extdata", "Amplitudes/", package="PoDCall")
## Read in data files
ampData <- PoDCall::importAmplitudeData(dataDirectory=path, skipLines=0)
## Calculate thresholds, metrics, concentrations
thresholdTable <- PoDCall::podcallThresholds(plateData=ampData)
print(thresholdTable)

#4.4 -podcallChannelPlot() 
## Read in data and threshold table
path <- system.file("extdata", "Amplitudes/", package="PoDCall")
ampData <- PoDCall::importAmplitudeData(dataDirectory=path, skipLines=0)
data("thrTable")
thresholdTable <- thrTable
## Select channel and well to plot
ch <- 1 # target channel
well_id <- names(ampData)[1] # First well in list
## Plot title
plotTitle <- paste0(well_id, ", Ch", ch)
## Create plot
PoDCall::podcallChannelPlot(channelData=ampData[[well_id]][,ch],
                   thr=thresholdTable[well_id, "thr_target"],
                   channel=ch,
                   plotId=plotTitle)

#4.5 podcallScatterplot()
## Read in data and threshold table
path <- system.file("extdata", "Amplitudes/", package="PoDCall")
ampData <- PoDCall::importAmplitudeData(dataDirectory=path, skipLines=0)
thresholdTable <- thrTable
## Select channel and well to plot
ch <- 1 # target channel
well_id <- names(ampData)[1] # First well in list
## Plot title
plotTitle <- paste0(well_id, ", Ch", ch)
## Create plot
PoDCall::podcallScatterplot(channelData=ampData[[well_id]][,ch],
                   thr=thresholdTable[well_id, "thr_target"],
                   channel=ch,
                   plotId=plotTitle)


# 4.6  podcallHistogram()
## Read in data and threshold table
path <- system.file("extdata", "Amplitudes/", package="PoDCall")
ampData <- PoDCall::importAmplitudeData(dataDirectory=path, skipLines=0)
thresholdTable <- thrTable
## Select channel and well to plot
ch <- 1 # target channel
well_id <- names(ampData)[1] # First well in list
## Plot title
plotTitle <- paste0(well_id, ", Ch", ch)
## Create plot
PoDCall::podcallHistogram(channelData=ampData[[well_id]][,ch],
                 thr=thresholdTable[well_id, "thr_target"],
                 channel=ch,
                 plotId=plotTitle)


#4.7  podcallMultiplot()
## Read in data and threshold table
path <- system.file("extdata", "Amplitudes/", package="PoDCall")
ampData <- PoDCall::importAmplitudeData(dataDirectory=path, skipLines=0)
thresholdTable <- thrTable
## Channel to plot
ch <- 1
## Create comparison plot
podcallMultiplot(plateData=ampData,
                 thresholds=thresholdTable[names(ampData),], 
                 channel=ch)


## Create comparison plot  but get the data used for making the plot
pDCdd <- podcallMultiplot_dd(plateData=ampData,
                 thresholds=thresholdTable[names(ampData),], 
                 channel=ch)

pDCdd[4]