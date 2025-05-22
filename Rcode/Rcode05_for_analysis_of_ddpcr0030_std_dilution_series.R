#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
#remove everything in the working environment, without a warning!!
rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
#define working directory
wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2022/ddPCR_qPCR_MST"
wd00 <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR"
#define input file  directory
wd01 <- "ddpcr_resultater"
#define input file  directory
wd01 <- "ddpcr_resultater"
setwd(wd00)
#define external input directory
wd01ext <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6/output05_stdcrv_plots_and_tables_from_Rcode"
wd01ext <- paste0(wd00,"/MONIS6_2021_data/output05_stdcrv_plots_and_tables_from_Rcode")
# paste together working directory and input file  directory
wd00_01 <- paste0(wd00,"/",wd01)
#grep for the '.csv' files in the directory, and place in a list
lst_fcsv <- list.files(wd00_01)[grep("\\.csv",list.files(wd00_01))]
#lst_fcsv <- lst_fcsv[grep("ddPCR0015",lst_fcsv)]
#lst_fcsv <- lst_fcsv[grep("ddpcr0017",lst_fcsv)]
lst_fcsv <- lst_fcsv[grep("ddpcr0030",lst_fcsv)]
inf02 <- lst_fcsv
# read in file from input file directory
df_ddP <- read.csv(paste0(wd00_01,"/",inf02))

# define file name for external input file with 
# qPCR standard curve efficiencies and slopes 
inf03 <- "out05_MONIS6_eDNA_std_crv_efficiencies.csv" 
wd03.1_inf03 <- paste(wd01ext,"/",inf03,sep="")
# read in file with qPCR standard curve efficiencies and slopes
#read in csv with qPCR data from MxPro
df_qef01 <- read.csv(wd03.1_inf03,header=T,sep=",")
# use dplyr to select columns and group by these columns, and 
# then take the max, the stddev, and stderr
df_qef02 <- df_qef01  %>%
  dplyr::select(spec.lat, intc2, lod.val,loq.val,rEffic,ampF ) %>% 
  dplyr::group_by(spec.lat) %>%
  #dplyr::summarise_each(funs(max, sd, se=sd(.)/sqrt(n())), intc2)
  #dplyr::summarise(max_lod=max(lod.val))
  dplyr::summarise_each(funs(min, max), c(lod.val,loq.val,rEffic))
# get the genus and the species names
genusNm <- sapply(strsplit(df_qef02$spec.lat," "), "[[", 1)
speciesNm <- sapply(strsplit(df_qef02$spec.lat," "), "[[", 2)
AbrNm <- paste0(substr(genusNm,1,3),substr(speciesNm,1,3))
# add back to data frame
df_qef02$AbrNm <- AbrNm
# define input file with plate setup in excel
# inf03 <- "plate_setup_ddpcr0017_Mnelei.xls"
# inf03 <- "plate_setup_ddpcr0015_Myaare.xls"
inf03 <- "plate_setup_ddpcr0030_stddiltest.xls"
inf03 <- paste0(wd00_01,"/",inf03)
# read in plate setup from excel flie , skip the first 25 rows
#pltstp <- readxl::read_xls(inf03,skip = 25)
pltstp <- readxl::read_xls(inf03,skip = 124)
# replace double underscores
pltstp$WellName <- gsub("__","_",pltstp$WellName)
# match to get sample name into ddPCR plate setup 
df_ddP$Sample.description.5 <- pltstp$WellName[match(df_ddP$Well,pltstp$WellNumber)]
df_ddP$Sample.description.6 <- df_ddP$Sample.description.5
df_ddP$Sample.description.7 <- df_ddP$Sample.description.5
df_ddP$Sample.description.6[grepl("_std_",df_ddP$Sample.description.6)] <- "STD"

df_ddP$Sample.description.5 <- sapply(strsplit(df_ddP$Sample.description.6,"_"), "[[", 3)
df_ddP$Sample.description.8 <- sapply(strsplit(df_ddP$Sample.description.6,"_"), "[[", 2)
df_ddP$Sample.description.9 <- sapply(strsplit(df_ddP$Sample.description.6,"_"), "[[", 1)
df_ddP$stdlvl <- df_ddP$Sample.description.5
df_ddP$stdlvl[is.na(df_ddP$stdlvl)] <- 0
df_ddP$stdlvl[!df_ddP$Sample.description.8=="STD"] <- NA
df_ddP$stdlvl <- as.numeric(df_ddP$stdlvl)
df_ddP$Sample.description.1 <- df_ddP$Sample.description.8
# View(df_ddP)
# spcAbbr <- df_ddP$Sample.description.7[df_ddP$Sample.description.8=="STD"]
# spcAbbr <- sapply(strsplit(spcAbbr,"_"), "[[", 1)
# df_ddP$spcAbbr <- df_ddP$Sample.description.6
df_ddP$spcAbbr <- df_ddP$Sample.description.9
#head(df_ddP,6)
# Read in file from MxPro qPCR analysis
#define working directory
wd00.1 ="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6"
wd00.1 = paste0(wd00,"/MONIS6_2021_data")
#wd00 <- "/Users/steenknudsen/Documents/Documents/NIVA_Ansaettelse_2020/NOVANA_proever_2018_2019"
#define input directory
wd01.1 <- "output02_merged_txtfiles_from_mxpro_for_MONIS6"
wd02.1 <- paste(wd00.1,"/",wd01.1,sep="")
#define inputfile
inpf01 <- "outfile02_merged_mxpro_csvfls_MONIS6.csv"
#paste path and input file together
if01wp <- paste(wd02.1,"/",inpf01,sep="")
#read in csv with qPCR data from MxPro
df_e01 <- read.csv(if01wp,header=T,sep=";")
#limit to only comprise one welltype
df_e02 <-  df_e01[df_e01$WellType=="Standard",]
# make the column numeric
df_e02$Quantitycopies <- as.numeric(df_e02$Quantitycopies)
# get the lod max value for each assay
df_e02$loq.val_max <- df_qef02$loq.val_max[match(df_e02$speciesabbr,df_qef02$AbrNm)]
# replace NAs with zeroes
df_e02$Quantitycopies[is.na(df_e02$Quantitycopies)] <- 0
# only include rows that have 'Quantitycopies' above the 'loq.val_max'
df_e02 <- df_e02[(df_e02$Quantitycopies > df_e02$loq.val_max),]
# load the dplyr library
require(dplyr)
#  I used 3 uL template, this means the original sample contains 1/3
df_e02$cppuL <- as.numeric(df_e02$Quantitycopies/3)
# head(df_e02,4)
# unique(df_e02$speciesabbr)
df_e02$cppuL[is.na(df_e02$cppuL)] <- 0
# make the column with 'Conc.copies.µL.' numeric
df_ddP$Conc.copies.µL. <- as.numeric(df_ddP$Conc.copies.µL.)
#make a column with log10 to copies plus 1
df_ddP$log10.Conc.copies.µL <- log10(df_ddP$Conc.copies.µL.+1)
df_ddP$l10cp <-  df_ddP$log10.Conc.copies.µL 
# calculate the total copy number for the total reaction volumne
df_ddP$totcpin25uLrxn <- df_ddP$Conc.copies.µL.*25
# I used 5 uL of template so the confidence intervals should be multiplied with 5
df_ddP$ccpuL <- df_ddP$Conc.copies.µL.*5
# replace any NAs with zeros for the copy number
df_ddP$ccpuL[is.na(df_ddP$ccpuL)] <- 0
# get unique species abbreviations
spNms <- unique(df_ddP$spcAbbr)
#define empty list
lstg <- list()
# set a running number
i <- 1
#iterate over unique species names
for (spcN in spNms)
{
  print(spcN)
  df_ddPsp <- df_ddP[df_ddP$spcAbbr==spcN,]
  # get copies for only std
  expstd <- df_ddPsp$stdlvl[df_ddPsp$Sample.description.8=="STD"]
  obsstd <- df_ddPsp$ccpuL[df_ddPsp$Sample.description.8=="STD"]
  # plot the two std dilutin series from ddPCR and qPCR
  #plot(log10(expstd),log10(obsstd))
  df_diff01 <- data.frame(x=obsstd,y=expstd)
  df_diff02 <- data.frame(x=obsstd,y=expstd)
  # https://stackoverflow.com/questions/9977686/how-to-remove-rows-with-any-zero-value
  ##Go through each row and determine if a value is zero
  row_sub = apply(df_diff02, 1, function(row) all(row !=0 ))
  ##Subset as usual
  df_diff02 <- df_diff02[row_sub,]
  # also check for values above 5e5
  row_sub = apply(df_diff02, 1, function(row) all(row <5e5 ))
  ##Subset as usual
  df_diff02 <- df_diff02[row_sub,]
  # also check for values above 5e5
  row_sub = apply(df_diff02, 1, function(row) all(row >=10 ))
  ##Subset as usual
  df_diff02 <- df_diff02[row_sub,]
  
  df_diff02.1 <- log10(df_diff02)
  
  #Replace NaN & Inf with NA
  df_diff02.1[is.na(df_diff02.1) | df_diff02.1=="Inf"] = NA
  # check the data frame you will use for making the linear model
  if (nrow(df_diff02.1)!=0)
  {
    # make linear model for the standard dilution series
    lmd2 <- lm(formula=y~x, data=df_diff02.1)
    # get the intercept of the linear model
    intcp2 <- lmd2$coefficients[1]
    # copy the intercept into a factor
    f2 <- 1/(10^intcp2)
    }
  # if no rows are available
    else 
    {
      f2 <- 0
    }
  # add to growing list
  lstg[[i]] <- f2[[1]]
  #increase th count by one
  i <- i+1
# end iteration over species
}
# unlist the factors
lstg1 <- unlist(lstg)
#bind the rows in each list in to one data frame
datbl_f2 <- cbind(lstg1,spNms, fill=T)
#make it a data fram
df_f2 <- as.data.frame(datbl_f2)
# copy column with species abbreviation name
df_f2$spNms2 <- df_f2$spNms  
# substitute in different species abbreavtion names
df_f2$spNms2 <- gsub("PsFa28" ,"Psefar" ,df_f2$spNms2)
df_f2$spNms2 <- gsub("PsVe28" ,"Psever" ,df_f2$spNms2)
spNms3 <- df_f2$spNms2
# unique(df_e02$speciesabbr)
# match back to get factor for conversion for qPCR measurements
df_e02$f2 <- df_f2$lstg1[match(df_e02$speciesabbr,df_f2$spNms2)]
df_e02$f2 <- as.numeric(df_e02$f2)
#View(df_e02)
# only retain standard dilution series
df_e02 <- df_e02[df_e02$WellType=="Standard",]
# qlod <- min(expstd[expstd>=1E1])
#factor
# modify the copy count in the qPCR 
df_e02$cppuLm <- df_e02$cppuL*df_e02$f2
# use dplyr to select columns and group by these columns, and
# then take the mean, the stddev, and stderr
df_e03 <- df_e02 %>%
  dplyr::select(Replicate, speciesabbr, cppuLm,smpltp,WellType ) %>% 
  dplyr::group_by(Replicate,speciesabbr,smpltp) %>%
  dplyr::summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), cppuLm)
# replace NA with zeroes in multiple columns
df_e03 <- df_e03 %>% 
  dplyr::mutate_at(c('mean','sd','se'), ~tidyr::replace_na(.,0))
# add a prefix to column names to indicate the data comes from qPCR results
colnames(df_e03) <- paste0("qP.",colnames(df_e03))

#make a column with log10 to copies plus 1
df_ddP$log10.Conc.copies.µL <- log10(df_ddP$Conc.copies.µL.+1)
df_ddP$l10cp <-  df_ddP$log10.Conc.copies.µL 
# calculate the total copy number for the total reaction volumne
#df_ddP$Conc.copies.µL.[df_ddP$stdlvl!=0]
df_ddP$totcpin25uLrxn <- df_ddP$Conc.copies.µL.*25
# get the log10 to the total count of molecules in the reaction 
#inferred by the ddPCR machine
df_ddP$l10trxn <- log10(df_ddP$totcpin25uLrxn+1)
# copy the column with sample number
df_ddP$smpl.No <- df_ddP$Sample.description.1
# copy the column with template volume
df_ddP$templvol <- df_ddP$Sample.description.2
# df_ddP$templvol <- as.numeric(gsub("voltempl","",df_ddP$Sample.description.2))
# df_ddP$templvol <- as.numeric(gsub("templvol","",df_ddP$Sample.description.2))
df_ddP$templvol <- as.numeric(gsub("templvol_([0-9]{+}).*","\\1",df_ddP$Sample.description.2))
# copy the column with standard dilution level
#df_ddP$stdlvl <- df_ddP$Sample.description.3
df_ddP$stdlvl[is.na(df_ddP$stdlvl)] <- 0

df_ddP$dsc2std <- df_ddP$stdlvl
#dsc2 <- as.character(log10(df_ddP$Sample.description.3))
#df_ddP$dsc2.2 <- gsub("^","STD5E",dsc2)
#df_ddP$smpl.No[df_ddP$smpl.No=="STD"] <- df_ddP$dsc2.2[df_ddP$smpl.No=="STD"]
# take log10 to the upper and lower confidence intervals of the concentration inferred
df_ddP$ddPBR_sdmn_l10 <- log10(df_ddP$PoissonConfMin+1)
df_ddP$ddPBR_sdmx_l10 <- log10(df_ddP$PoissonConfMax+1)
# take log10 to the upper and lower confidence intervals of the 
#total copy number in the reaction
df_ddP$l10PcMi_ttcnt <- log10((df_ddP$PoissonConfMin*25)+1)
df_ddP$l10PcMx_ttcnt <- log10((df_ddP$PoissonConfMax*25)+1)
# assign column for volume of template added
#df_ddP$tvadd <- as.numeric(gsub("templvol(.*)uL","\\1",df_ddP$Sample.description.3))
df_ddP$tvadd <- df_ddP$templvol
# take log10 to the expected concentration, and divide by the total reaction volume
df_ddP$l10.dsc2 <- log10(df_ddP$dsc2std/25)
# take log10 to the expected copy number , DO NOT  divide by the total reaction volume
df_ddP$l10.tcrx_added <- log10(df_ddP$dsc2std)
# replace the NAs with 0
df_ddP$l10.dsc2[is.na(df_ddP$l10.dsc2)] <- 0
df_ddP$l10.tcrx_added[is.na(df_ddP$l10.tcrx_added)] <- 0

df_ddP$l10.dsc2[is.infinite(df_ddP$l10.dsc2)] <- 0
df_ddP$l10.tcrx_added[is.infinite(df_ddP$l10.tcrx_added)] <- 0
# min(df_ddP$Accepted.Droplets)
# max(df_ddP$Accepted.Droplets)
#_________________________________________________________________________________

# copy data frame
df_e04 <- df_e03
#View(df_e03)
# grep for the standard dilution
dPCRstdmspls <- df_ddP$smpl.No[grepl("STD",df_ddP$smpl.No)]
# replace capital letters with small letters
dPCRstdmspls <- tolower(dPCRstdmspls)
# substitute to retain standard dilution factor but not volume
dPCRstdmspls <-gsub("std(.)e(.*)","stdE\\2",dPCRstdmspls)
# use to replace in the original data frame
df_ddP$smpl.No[grepl("STD",df_ddP$smpl.No)] <- dPCRstdmspls
# grep for the standard dilution samples in the qPCR data
qPCRstdmspls <- df_e04$qP.smpltp[grepl("std",df_e04$qP.smpltp)]
# substitute to retain standard dilution factor but not volume
qPCRstdmspls <-gsub("std(.)E(.*)","stdE\\2",qPCRstdmspls)
# grep for the 'NTC' sample to replace this sample name
qPCRstdmspls[grepl("NTC",qPCRstdmspls)] <- "NTC"
# add back the substituted qPCR standard dilution names
df_e04$qP.smpltp[grepl("std",df_e04$qP.smpltp)] <- qPCRstdmspls
# Check only the standard dilution positive controls in the ddPCR setup
dfg <- df_ddP[df_ddP$Sample.description.1=="STD",]
# substitute in the std name and paste together with lowercase letters
dfg$Sample.description.1 <- paste0(tolower(dfg$Sample.description.1), gsub("^(.)E(.*)$","E\\2",dfg$Sample.description.5))
# replace the column with smpl.No for the rows in 'Sample.description.1 that match 'STD'   
df_ddP$smpl.No[df_ddP$Sample.description.1=="STD"] <- dfg$Sample.description.1

# I used 5 uL of template so the confidence intervals should be multiplied with 5
df_ddP$ccpuL <- df_ddP$Conc.copies.µL.*5
# find sampletypes included in ddPCR setup not included in qPCR setup
# and addthem to vector
extsmplTp <- df_ddP$smpl.No[!df_ddP$smpl.No %in% df_e04$qP.smpltp]
# get the 'tidyverse' functions
library(tidyverse)
# add rows to the data frame: https://stackoverflow.com/questions/28467068/how-to-add-a-row-to-a-data-frame-in-r
df_e04.1 <- as.data.frame(df_e04) %>% add_row("qP.smpltp" = extsmplTp)
# make it a tibble again
df_e04 <- as_tibble(df_e04.1)
# subset to only include the species tested in the ddPCR
df_e04 <- df_e04[df_e04$qP.speciesabbr %in% spNms3 ,]
# exclude the test for Rhihar as this failed in this setup
#df_e04 <- df_e04[df_e04$qP.speciesabbr="Rhihar" ,]
#View(df_e04)

df_ddP$spcAbbr2 <- df_f2$spNms2[match(df_ddP$spcAbbr,df_f2$spNms)]
df_ddP$sA_sN <- paste0(df_ddP$spcAbbr2,"_",df_ddP$smpl.No)
df_e04$sA_sN <- paste0(df_e04$qP.speciesabbr,"_",df_e04$qP.smpltp)
# split text string by using underscore as delimiter, and get third element
stdlvl2 <- sapply(strsplit(df_ddP$Sample.description.7,"_"), "[[", 3)
# get the species abbreviations
spcAbbr3 <- df_ddP$spcAbbr

# bind columns and make it a tibble
# use this tibble to later on replace the 'df_e04' data frame
# that you use for making the ggplot 
tibl_ddP2 <- as_tibble(cbind(spcAbbr3,stdlvl2))
# change the column names
colnames(tibl_ddP2) <- c("dP.speciesabbr", "dP.smpltp")
# add copy count per uL
tibl_ddP2$dP.ccpuL <- df_ddP$ccpuL
#replace any NAs with zeroes
tibl_ddP2$dP.ccpuL[is.na(tibl_ddP2$dP.ccpuL)] <- 0
# modify the 'PoissonConfidence' intervals
# I used 5 uL of template so the confidence intervals should be multiplied with 5
df_ddP$PCmn68t5 <- df_ddP$PoissonConfidenceMin68*5
df_ddP$PCmx68t5 <- df_ddP$PoissonConfidenceMax68*5
# and then add them to the tibble
tibl_ddP2$dP.pcmn <- df_ddP$PCmn68t5
tibl_ddP2$dP.pcmx <- df_ddP$PCmx68t5
#replace any NAs with zeroes
tibl_ddP2$dP.pcmn[is.na(tibl_ddP2$dP.pcmn)] <- 0
tibl_ddP2$dP.pcmx[is.na(tibl_ddP2$dP.pcmx)] <- 0

# match between data frame to combine qPCR and ddPCR results
df_e04$dP.ccp <- df_ddP$ccpuL[match(df_e04$sA_sN,df_ddP$sA_sN)]

# copy the standard deviation for the qPCR results to have a minimum and a maximum
# to match the ddPCR 
df_e04$qP.sdmn <- df_e04$qP.mean-df_e04$qP.sd
df_e04$qP.sdmx <- df_e04$qP.mean+df_e04$qP.sd
# match to get confidence intervals for ddPCR evaluations
df_e04$dP.pcmn <- df_ddP$PoissonConfMin[match(df_e04$sA_sN,df_ddP$sA_sN)]
df_e04$dP.pcmx <- df_ddP$PoissonConfMax[match(df_e04$sA_sN,df_ddP$sA_sN)]

# I used 5 uL of template so the confidence intervals should be multiplied with 5
df_ddP$PCmn68t5 <- df_ddP$PoissonConfidenceMin68*5
df_ddP$PCmx68t5 <- df_ddP$PoissonConfidenceMax68*5
# match back to get upper and lower confidence ingterval
df_e04$dP.pcmn <- df_ddP$PCmn68t5[match(df_e04$sA_sN,df_ddP$sA_sN)]
df_e04$dP.pcmx <- df_ddP$PCmx68t5[match(df_e04$sA_sN,df_ddP$sA_sN)]

# replace NA with zeroes in multiple columns
df_e04 <- df_e04 %>% 
  mutate_at(c('dP.ccp','dP.pcmn','dP.pcmx'), ~replace_na(.,0))
# 
# df_e04$dP.pcmn <- log10(df_e04$dP.pcmn+1)
# df_e04$dP.pcmx <- log10(df_e04$dP.pcmx+1)
# df_e04$dP.ccp <- log10(df_e04$dP.ccp+1)
# df_e04$qP.sdmn <- log10(df_e04$qP.sdmn+1)
# df_e04$qP.sdmx <- log10(df_e04$qP.sdmx+1)
# df_e04$qP.mean <- log10(df_e04$qP.mean+1)

# use select in dplyr to keep specified columns. dplyr requires a grouping variable
df_e04 <- df_e04 %>%  dplyr::group_by(qP.smpltp,sA_sN) %>% 
  dplyr::select(sA_sN,qP.smpltp,qP.mean,dP.ccp,qP.sdmn,qP.sdmx,dP.pcmn,dP.pcmx) 
# rename column headers
colnames(df_e04) <- c("spcAbbr","smplNm","qP.mcp",
                      "dP.mcp","qP.sdmn","qP.sdmx","dP.sdmn","dP.sdmx")
# get only first part of string split using underscore as delimeter
df_e04$spcAbbr <- sapply(strsplit(df_e04$spcAbbr,"_"), "[[", 1)
#View(df_e04)
# use gather in tidyr package to  
df_e05 <- df_e04 %>% tidyr::gather(key = "cattype", value = "cpcnt", -c(spcAbbr ,smplNm))
# use separate to have two colunms, one for the machne type, and one for the category of
# counts
df_e05 <- tidyr::separate(data = df_e05, col = cattype, 
                          into = c("machine", "categ"), sep = "\\.")
#Remove rows with NA
df_e05.1  <- df_e05[(df_e05$spcAbbr!="NA"),]
# use 'pivot_wider' from the 'tidyr' package
# like here :https://stackoverflow.com/questions/30592094/r-spreading-multiple-columns-with-tidyr
df_e05.1 <- tidyr::pivot_wider(data=df_e05.1,id_cols = c(spcAbbr,smplNm, machine), 
                   names_from = categ, 
                   values_from = cpcnt)
# copy tibble back
df_e05 <- df_e05.1
# substitute to replace abbreviated machine names with long machine names
df_e05$machine <- gsub("dP","ddPCR",df_e05$machine)
df_e05$machine <- gsub("qP","qPCR",df_e05$machine)
# make a list with colors for the sampling categories
clfH <- rep(c("grey54","white"),length(unique(df_e05$smplNm)))
# count the number of elements
ncl <- length(clfH)
clfH[(ncl-7):(ncl)] <- rep(c("orchid4","orchid1"),4)
# use only purple colors for the standard dilution series
clfH <- rep(c("orchid4","orchid1"),9)
# define a list with sample ids for control samples, in the order we want them to appear in the figure
controlnames <- c("NTC",
                  "stdE-2",
                  "stdE-1",
                  "stdE0",
                  "stdE1",
                  "stdE2",
                  "stdE3",
                  "stdE4",
                  "stdE5")
# get unique sample ids from the dataframe used for plotting
listSmplNm <- unique(df_e05$smplNm)
# exclude control names from the list of unique sample ids
listSmplNm <-    listSmplNm[!listSmplNm %in% controlnames]
# add control names back to the list of unique sample ids (at the end), now in the order we want
listSmplNm <-    c(listSmplNm, controlnames) 
# convert smplNm to a factor with the levels (order) given by the list we have made
df_e05$smplNm <- factor(df_e05$smplNm, levels= listSmplNm)
# remember that for ggplot() we should now use “fill=smplNm” instead of “fill=as.factor(smplNm)”

# add a column for the machine type
tibl_ddP2$machine <- "ddPCR"
# change the column names
colnames(tibl_ddP2) <- c("spcAbbr","smplNm","mcp","sdmn","sdmx","machine")

tibl_ddP2 <- tibl_ddP2[!tibl_ddP2$mcp>2E5,]
tibl_ddP2 <- tibl_ddP2[!tibl_ddP2$mcp<2E-1,]
# copy the tibble
df_e06 <- tibl_ddP2
#View(df_e05)
#df_e05[df_e05$spcAbbr=="Psever",]
# https://stackoverflow.com/questions/41047939/ggplot-barplot-how-to-display-small-positive-numbers-with-log-scaled-y-axis
pd = position=position_dodge(.5)
plt07 <- ggplot(data=df_e06, aes(x=smplNm,y=mcp)) +
  theme_bw() +
  theme_classic() +
  facet_wrap(~spcAbbr, ncol =4) +
  # Der var 2 ting der skulle fikses og så en tredje ting jeg tror er en fejl
  # 1) y-position af tiles
  # geom_tile skal bruge x og y ”mappings” eller ”aesthetics” (aes). Som du nok ved, fordi de er ikke specificerede i geom_tile(), bliver der anvendt dem fra ggplot() kaldet.
  # Dvs. x=as.factor(smplNm), y=mcp
  # Problemet er at nogle af mcp værdier er lige med 0 og derfor ligger centerpunkterne for de tilhørende tiles uden for din figur. Man kan diskutere om en tile med uendelig højde burde alligevel strække tilbage og ind over figuren men den må vi tage en anden dag. Den nemme løsning var at sætte y=1, eller en anden værdi som ligger indenfor den synlige y-range af din figur.
  # 2) farveskema – dette havde du sikkert selv fundet ud af efter at 1) var fikset.
  # geom_tile() bruger du ”fill=clfH”. Her har du brugte paletten. ”fill” skal referere til den variabel der bestemmer farven og ikke paletten der skal anvendes. Derefter mangler du scale_fill_manual. Så nu ser geom_tile sådan her ud:
  # geom_tile(aes(height=Inf,fill=as.factor(smplNm),y=1,alpha=0.4),
  #           position=position_nudge(x=-0.05), show.legend=F) +
  # scale_fill_manual(values=clfH) +

geom_tile(aes(height=Inf,fill=as.factor(smplNm),y=1), 
          position=position_nudge(x=-0.05), show.legend=F) +
  scale_fill_manual(values=alpha(clfH,0.04)) +
  
  scale_y_log10(limits=c(1e-2, 1e5),
                breaks = scales::trans_breaks("log10", function(y) 10^y),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(panel.background = element_rect(fill = "transparent")) +
  # add points and error bars
  geom_errorbar(aes(ymin=sdmn, ymax=sdmx, colour=machine), width=0.9, lty=1,position=pd, size=0.8) +
  geom_point(aes(colour=machine), position=pd, size=3) +
  # add horizontal lines to show decades on y-axis
  geom_hline(yintercept=1E5, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E4, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E3, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E2, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E1, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E0, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E-1, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E-2, color=scales::alpha("gray53",0.4), size=0.12) +
  # add line for LOQ and/ or LOD
  #geom_hline(yintercept=2E-1, color="springgreen4") +
  #add a text label to the horizontal line
  # geom_text(aes(x = (as.factor(smplNm)[90]), y = 2E0-0.5E0), size=2.4,
  #           label = "nedre grænse for kvantificering", color = "springgreen4") +
  
  #geom_hline(yintercept=3*10^1, color="#E69F00",lty=2) +
  #geom_hline(yintercept=qlod*f2, color="#E69F00",lty=1) +
  theme(axis.text.x = element_text(angle = 90, hjust=1))   +
  ylab("målt konc eDNA molekyler/ uL ") + xlab("forventet konc") +
  annotation_logticks(sides="l") +
  #coord_flip()+
  scale_color_manual(values=c("springgreen4", "#E69F00"))
plt07  <- plt07 + labs(title = "standard fortyndings række i ddPCR og qPCR")#,

#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
plt07 <- plt07 + labs(color='maskine')
plt07 <- plt07 + labs(fill='maskine')
plt07 <- plt07 + labs(shape='maskine')

plt07t <- plt07 + labs(title = " ")#,
plt07t
# modify the input file name
inf02_1 <- gsub("\\.csv","",inf02)
#make filename to save plot to save the plot to
fgNm07 <- paste0(wd00_01,"/Fig02_v01",inf02_1,".png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- TRUE
if(bSaveFigures==T){
  ggsave(plt07t,file=fgNm07,#width=210,height=297,
         width=297,height=210,
         units="mm",dpi=600)
}

# colnames(df_e05)
#_______________________________________________________________________________
# example with tidyr
#_______________________________________________________________________________
library(tidyr)
# L1 <- rep(LETTERS[1:4],2)
# L2 <- c(rep(1,4),rep(2,4))
# N1 <- round(runif(8),2)
# N2 <- round(runif(8),2)
# N3 <- round(runif(8),2)
# df01 <- as.data.frame(cbind(L1,L2,N1,N2,N3))
# df01 <- df01[order(df01$L1),]
# df01 %>%
#   tidyr::pivot_wider(names_from="L2",
#   values_from=c("N1","N2","N3"))
# #_______________________________________________________________________________
# use tidyr example above
df_e06  <- df_e05 %>%
  tidyr::pivot_wider(names_from="machine",
                     values_from=c("mcp","sdmn","sdmx"))
library(ggplot2)
df_e06$smplTp <- "MST"
df_e06$smplTp[grepl("std",df_e06$smplNm)] <- "STD"
plt08 <- ggplot(data=df_e06, aes(x=mcp_qPCR ,y=mcp_ddPCR) ) +
  theme_classic() +
  facet_wrap(~spcAbbr, ncol =4) +
  geom_point(aes(fill=smplTp, shape=smplTp), size=3) + 
  scale_shape_manual(values = c(21)) +
  scale_fill_manual(values=c("pink")) +
  geom_errorbar(aes(ymin =  sdmn_ddPCR,ymax =  sdmx_ddPCR)) + 
  geom_errorbarh(aes(xmin =sdmn_qPCR,xmax = sdmx_qPCR)) +
  scale_y_log10(limits=c(1e-2, 1e5),
                breaks = scales::trans_breaks("log10", function(y) 10^y),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(panel.background = element_rect(fill = "transparent")) +
  scale_x_log10(limits=c(1e-2, 1e5),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(panel.background = element_rect(fill = "transparent")) +
  # add horizontal lines to show decades on y-axis
  geom_hline(yintercept=1E5, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E4, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E3, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E2, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E1, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E0, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E-1, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_hline(yintercept=1E-2, color=scales::alpha("gray53",0.4), size=0.12) +
  # add horizontal lines to show decades on y-axis
  geom_vline(xintercept=1E5, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_vline(xintercept=1E4, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_vline(xintercept=1E3, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_vline(xintercept=1E2, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_vline(xintercept=1E1, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_vline(xintercept=1E0, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_vline(xintercept=1E-1, color=scales::alpha("gray53",0.4), size=0.12) +
  geom_vline(xintercept=1E-2, color=scales::alpha("gray53",0.4), size=0.12) +
  xlab("eDNA molekyler /uL i qPCR") + ylab("eDNA molekyler /uL i ddPCR") +
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, color="blue") +
  labs(color='type') +
  labs(fill='type') +
  labs(shape='type') +
  #labs(title = "sandmusling - 'Mya arenaria' fra 2021")#,
  #labs(title = "ribbegoble - 'Mnemiopsis leidyi' fra 2021")#,
labs(title = "sammenligning mellem ddPCR og qPCR")#

plt08
fgNm08 <- paste0(wd00_01,"/Fig02_v02",inf02_1,".png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- TRUE
if(bSaveFigures==T){
  ggsave(plt08,file=fgNm08,#width=210,height=297,
         width=297,height=210,
         units="mm",dpi=600)
}

mchty <- c("ddPCR","qPCR") 
MSTsmpl071 <- df_e05$mcp[grepl("MST2021-071",df_e05$smplNm)]
MSTsmpl058 <- df_e05$mcp[grepl("MST2021-058",df_e05$smplNm)]
MSTsmpl085 <- df_e05$mcp[grepl("MST2021-085",df_e05$smplNm)]
df_MSTres <- as.data.frame(rbind(mchty,MSTsmpl071,MSTsmpl058,MSTsmpl085))
colnames(df_MSTres) <- df_MSTres[1,]
df_MSTres <- df_MSTres[-1,]

# ------------- plot Combined figure -------------
# USe patchwork to make combined figure
library(patchwork)
# set a variable to TRUE to determine whether to save figures
bSaveFigures <- T
#getwd()
#define a filename to save to
fnm02 <- paste0(wd00_01,"/Fig0030_",inf02_1,".png")
#see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
# on how to arrange plots in patchwork

plt07t <- plt07 + labs(title = "A")#,
plt08t <- plt08 + labs(title = "B")#,
p <-  plt07t  +
  plt08t +
  
  plot_layout(nrow=1,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  patchwork::plot_annotation(caption="Fig0030_v01_stddilcomp") #& theme(legend.position = "bottom")
#p
if(bSaveFigures==T){
  ggsave(p,file=fnm02,width=297*1.2,height=210*0.8,
         units="mm",dpi=300)
}