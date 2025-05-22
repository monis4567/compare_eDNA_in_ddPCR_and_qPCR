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
wd00  <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR"
#define input file  directory
wd01 <- "ddpcr_resultater"
#define input file  directory
wd01 <- "ddpcr_resultater"
setwd(wd00)
getwd()
# paste together working directory and input file  directory
wd00_01 <- paste0(wd00,"/",wd01)
#grep for the '.csv' files in the directory, and place in a list
lst_fcsv <- list.files(wd00_01)[grep("\\.csv",list.files(wd00_01))]
# grep only files that begin with ddpcr and 4 digits 
# and 1 capital letter and 5 small letters 
lst_fcsv <- lst_fcsv[grepl("^ddpcr[0-9]{4}_[A-Z]{1}[a-z]{5}",lst_fcsv, ignore.case = F)]
# add wd to all files
lst_inf02 <- wd00_01 %>% paste0("/",lst_fcsv)
# https://statisticsglobe.com/merge-csv-files-in-r
library(readr)
# To read in all files from the list, and also make sure all columns are 
# characters, and the bind rows in to a tibble
tibl_csvs01 <-
  lst_inf02 %>%
  lapply(read_csv, col_types=cols(.default = col_character())) %>%
  bind_rows
# copy the sample description columns
tibl_csvs01$smpldscr01 <- tibl_csvs01$`Sample description 1`
tibl_csvs01$smpldscr02 <- tibl_csvs01$`Sample description 2`
tibl_csvs01$smpldscr03 <- tibl_csvs01$`Sample description 3`
tibl_csvs01$smpldscr04 <- tibl_csvs01$`Sample description 4`
# get the rows in the 'Target' column, that does not have the 6 letter 
# species abbreviation
tibl_csvs01$Target[!grepl("[A-z]",tibl_csvs01$Target)] <- tibl_csvs01$smpldscr01[!grepl("[A-z]",tibl_csvs01$Target)]
# get the 6 letter name for the assay applied , loose any extra characters beoynd 6
tibl_csvs01$smpldscr04<- substring(tibl_csvs01$Target,1,6)
# make the 'MST' samples in the 'smpldscr01' column appear as "unknown"
tibl_csvs01$smpldscr01[grepl("MST",tibl_csvs01$smpldscr01)] <- "unknown"
tibl_csvs01$smpldscr01[grepl("NEKFeb",tibl_csvs01$smpldscr01)] <- "unknown"
# change for the 6 letter sample abbreviation
tibl_csvs01$smpldscr01[grepl("^[A-z]{6}$",tibl_csvs01$smpldscr01)] <- tibl_csvs01$smpldscr03[grepl("^[A-z]{6}$",tibl_csvs01$smpldscr01)]
# replace the NAs
tibl_csvs01$smpldscr01[is.na(tibl_csvs01$smpldscr01)] <- "unknown"
# edit the column for the standard dilution level
tibl_csvs01$smpldscr03 <- NA
tibl_csvs01$smpldscr03[grepl("A01",tibl_csvs01$Well)] <- "1E5"
tibl_csvs01$smpldscr03[grepl("B01",tibl_csvs01$Well)] <- "1E4"
tibl_csvs01$smpldscr03[grepl("C01",tibl_csvs01$Well)] <- "1E3"
tibl_csvs01$smpldscr03[grepl("D01",tibl_csvs01$Well)] <- "1E2"
tibl_csvs01$smpldscr03[grepl("E01",tibl_csvs01$Well)] <- "1E1"
tibl_csvs01$smpldscr03[grepl("F01",tibl_csvs01$Well)] <- "1E0"
tibl_csvs01$smpldscr03[grepl("G01",tibl_csvs01$Well)] <- "1E-1"
tibl_csvs01$smpldscr03[grepl("H01",tibl_csvs01$Well)] <- "0"
# add a column for vol of template added
tibl_csvs01$smpldscr02 <- "voltempl_5uL"
# copy the tibble as a data frame
df_ddP  <- as.data.frame(tibl_csvs01)
# define input file with plate setup in excel
inf03 <- "plate_setup_ddpcr0023_Psefar.xls"
# 
inf03 <- paste0(wd00_01,"/",inf03)
# read in plate setup from excel flie , skip the first 25 rows
pltstp <- readxl::read_xls(inf03,skip = 25)
# grep in the columns that '_std_' in the name, and replace in these names 
pltstp$WellName[grepl("_std_",pltstp$WellName)] <- 'STD'
# match to get sample name into ddPCR plate setup 
df_ddP$Sample.description.5 <- pltstp$WellName[match(df_ddP$Well,pltstp$WellNumber)]
df_ddP$Sample.description.6 <- df_ddP$Sample.description.5

df_ddP$Sample.description.6[grepl("_std_",df_ddP$Sample.description.6)] <- "STD"

stds <- df_ddP$Sample.description.5[df_ddP$smpldscr01=="STD"]
# add the dilution level of standard to a column 
df_ddP$stdlvl <- df_ddP$smpldscr03
df_ddP$stdlvl[is.na(df_ddP$stdlvl)] <- 0
df_ddP$stdlvl[!df_ddP$Sample.description.6=="STD"] <- NA
df_ddP$stdlvl <- as.numeric(df_ddP$stdlvl)
# multiply by 5 because you added 5 uL of the std dilution series
#df_ddP$stdlvl <- df_ddP$stdlvl*5
df_ddP$Sample.description.1 <- df_ddP$Sample.description.6
# Read in file from MxPro qPCR analysis
#define working directory
#wd00.1 ="/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6"
wd00.1 =paste0(wd00,"/MONIS6_2021_data")
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

# copy data frame
df_e02 <-  df_e01
# make the column numeric
df_e02$Quantitycopies <- as.numeric(df_e02$Quantitycopies)
#_______________________________________________________________________________
#___________________________________
# calculate LOD and LOQ for qPCR data
#___________________________________
df_qPCRm01 <- df_e02
df_qPCRm01$snm_pltn <- paste0(df_qPCRm01$speciesabbr,"_",df_qPCRm01$plateno)
df_qPCRm01$qcp <- df_qPCRm01$Quantitycopies
df_qPCRm01$qcp[is.na(df_qPCRm01$qcp)] <- 0
df_qPCRm01 <- df_qPCRm01[df_qPCRm01$qcp>0,]
#get the LOD for each qPCR run
#use the function aggregate to get the minimum value for a group
lodtable1 <- aggregate(df_qPCRm01[, "qcp"], list(df_qPCRm01$snm_pltn, df_qPCRm01$WellType), min)
#subset this table by group
lodtable2 <- lodtable1[ which(grepl("Standard",lodtable1$Group.2)), ] # use this line for the MxPro Machine
#rename the column names - last column holds 'Limit Of Detection' per
# qPCR run
colnames(lodtable2) <- c("spc.pltno","WellT","LOD")
#subset the data frame to only incl standard dilutions
# this time from the data frame that includes the 'failed' 
# standard - i.e the standards that di not amplify
fg <- df_qPCRm01[ which(df_qPCRm01$smpltp=="Std"), ]
fg <- df_qPCRm01[ which(df_qPCRm01$WellType=="Standard"), ]
#count the number of replicate used per plate
#see this webpage: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
#and this webpage: https://stackoverflow.com/questions/9809166/count-number-of-rows-within-each-group
lodtb3 <- fg %>% dplyr::count(snm_pltn,smpltp)
# get number of intended replicates
lodtb3$nint.rpl<- max(lodtb3$n)
#lodtb3[grepl("Neogo",lodtb3$gen_specnm.pltn),]
# paste together  species name and plate number and welltype content
lodtable2$nrepl <- lodtb3$n[match(lodtable2$spc.pltno, lodtb3$snm_pltn)]
#Now identify LOQ for each qPCR run
#limit the dataframe to only well type that equals standard
oc <- df_qPCRm01[(df_qPCRm01$WellType=='Standard'),] #use with MxPro
#exclude samples in the 'Std' Standard curve that did no amplify
oc <- oc[!is.na(as.numeric(as.character(oc$CtdRn))),] #use with MxPro
oc <- oc[oc$CtdRn!=0,]
#add a new column that merges two columns for species and qPCR plate no
oc$Quancp.spc.pltn <- paste(oc$Quantitycopies, oc$snm_pltn,  sep=".")
#count the occurences of dilution steps - i.e. the number of succesful replicates
#see this webpage: https://www.miskatonic.org/2012/09/24/counting-and-aggregating-r/
#and this webpage: https://stackoverflow.com/questions/9809166/count-number-of-rows-within-each-group
od <- oc %>% dplyr::count(Quantitycopies, Quancp.spc.pltn)
#turn this into a dataframe
oe<-as.data.frame(od)
#The 'n' column now holds the count of wells at this specific
# dilution level that successfully amplified
# lowest number below the number of total replicates will be
# the level just below the LOQ
#add a new column that merges two columns
#match the dilution step to the number of occurences -i.e. match between the two dataframes
no.occ <- oe$n[match(oc$Quancp.spc.pltn,oe$Quancp.spc.pltn)]
#add this column with counted occurences to the limited dataframe
og <- cbind.data.frame(oc,no.occ)
#get the number of replicates used
og$nrepl <- lodtable2$nrepl[match(og$snm_pltn,lodtable2$spc.pltno)]
#exlude all observations where 
#less than '3'number of replicates' amplified
oh<-og[(og$no.occ>=og$nrepl),]
oh$Quantitycopies <- as.numeric(oh$Quantitycopies)
#get the lowest dilution step that succesfully amplified on all 3 repliactes
#use aggregate to get the minimum for each
loqtable1 <- aggregate(oh[, "Quantitycopies"], list(oh$snm_pltn), min)
#change the column names
colnames(loqtable1) <- c("spc.plt","LOQ")
#copy the LOD table and add the corresponding LOQ values
loq.lod.table <- lodtable2
loq.lod.table$LOQ <- loqtable1$LOQ[match(lodtable2$spc.pltno,loqtable1$spc.plt)]
#append limit of quantification back to main data frame
df_qPCRm01$LOQ <- loq.lod.table$LOQ[match(df_qPCRm01$snm_pltn,loq.lod.table$spc.pltno)]
#append limit of detection back to main data frame
df_qPCRm01$LOD <- loq.lod.table$LOD[match(df_qPCRm01$snm_pltn,loq.lod.table$spc.pltno  ) ]
df_qPCRm01$LOQ <- as.numeric(df_qPCRm01$LOQ)
df_qPCRm01$LOD<- as.numeric(df_qPCRm01$LOD)
#make the 'Quantitycopies' numeric, and replace NAs with zeroes
df_qPCRm01$QuanCp02 <- as.numeric(df_qPCRm01$Quantitycopies)
df_qPCRm01$QuanCp02[is.na(df_qPCRm01$QuanCp02)] <- 0

loq.lod.table$spcabbr <- sapply(strsplit(loq.lod.table$spc.pltno,"_"), "[[", 1)
# get mean for the LOD and LOQ per species in qPCR
df_qP.lod.loq <- loq.lod.table %>%
  dplyr::select(spcabbr, LOD,LOQ ) %>% 
  dplyr::group_by(spcabbr) %>%
  dplyr::summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), c(LOD,LOQ))

#count by speciesabbr, smpltp, plateno
df_qPCRm02 <- df_qPCRm01 %>% 
  dplyr::group_by(speciesabbr, smpltp, plateno) %>% 
  dplyr::summarise(n = n()) 
# reduce to only get 'std' samples
df_qPCRm02.std <- df_qPCRm02[grepl("std",df_qPCRm02$smpltp),]
# get std dilution level as numeric value
df_qPCRm02.std$smpltp.dil <- as.numeric(gsub("std","",df_qPCRm02.std$smpltp))
# make the NAs 0
df_qPCRm02.std$smpltp.dil[is.na(df_qPCRm02.std$smpltp.dil)] <- 0
# get number of maximum replicates
df_qPCRm02.mxrpl <- df_qPCRm02.std %>%
  dplyr::select(speciesabbr, smpltp,plateno,n ) %>% 
  dplyr::group_by(speciesabbr,plateno) %>%
  dplyr::summarise_each(max, c(n))
# paste to have common speceies abbrevitation and plate number
df_qPCRm02.std$spAb_plN <- paste0(df_qPCRm02.std$speciesabbr,"_",df_qPCRm02.std$plateno)
df_qPCRm02.mxrpl$spAb_plN <- paste0(df_qPCRm02.mxrpl$speciesabbr,"_",df_qPCRm02.mxrpl$plateno)
# match to get max number of replicates
df_qPCRm02.std$mxrpl <- df_qPCRm02.mxrpl$n[match(df_qPCRm02.std$spAb_plN,df_qPCRm02.mxrpl$spAb_plN)]
# exclude replicate sets with too low number of replicates 
df_qPCRm03 <- df_qPCRm02.std[(!df_qPCRm02.std$n<df_qPCRm02.std$mxrpl),]
# get minimum dilution level that amplified on 
df_qPCRm04 <- df_qPCRm03 %>%
  dplyr::select(speciesabbr, smpltp,plateno,n ,smpltp.dil) %>% 
  dplyr::group_by(speciesabbr) %>%
  dplyr::summarise_each(min, c(smpltp.dil))
#
#View(df_qPCRm02.std)

#_______________________________________________________________________________

# load the dplyr library
require(dplyr)
#  I used 3 uL template, this means the original sample contains 1/3
df_e02$cppuL <- as.numeric(df_e02$Quantitycopies/3)
# only retain detections in qPCR below 42 CtdRn
df_e02$CtdRn[df_e02$CtdRn=="NoCt"] <- 0
df_e02$CtdRn <- as.numeric(df_e02$CtdRn)
df_e02 <- df_e02[df_e02$CtdRn<42,]
#replace NAs with 0
df_e02$cppuL[is.na(df_e02$cppuL)] <- 0
# exclude negative controls above zero
# https://sparkbyexamples.com/r-programming/r-subset-data-frame-by-column-value/
# Subset Rows by Checking values on Multiple Columns
df_e02.1 <- df_e02[!((grepl('NTC',df_e02$smpltp)) & df_e02$cppuL > 0),]
df_e02.1 <- df_e02.1[!((grepl('NEK',df_e02.1$smpltp)) & df_e02.1$cppuL > 0),]
# replace old version of  data frame
df_e02 <- df_e02.1
# copy column
df_ddP$Conc.copies.µL. <- df_ddP$`Conc(copies/µL)`
#make copy count
df_ddP$Conc.copies.µL. <- as.numeric(df_ddP$Conc.copies.µL.)
# replace NAs with zeroes
df_ddP$Conc.copies.µL.[is.na(df_ddP$Conc.copies.µL.)] <- 0
#make a column with log10 to copies plus 1
df_ddP$log10.Conc.copies.µL <- log10(df_ddP$Conc.copies.µL.+1)
df_ddP$l10cp <-  df_ddP$log10.Conc.copies.µL 
# calculate the total copy number for the total reaction volumne
#df_ddP$Conc.copies.µL.[df_ddP$stdlvl!=0]
df_ddP$totcpin25uLrxn <- df_ddP$Conc.copies.µL.*25
# I used 5 uL of template so the confidence intervals should be multiplied with 5
df_ddP$ccpuL <- df_ddP$Conc.copies.µL.*5
# make empty lists to add to
lst_f2 <- list()
lst_ddlod <- list()
lst_intcp2 <- list()
# iterate over species
lst_spcAb <- unique(df_ddP$smpldscr04)
# make a number to iterate over
i <- 1
#match("Karmik",lst_spcAb)
# iterate over species
for (sp in lst_spcAb)
{
  # see the species name
  print(sp)
  #}
  # subset the data frame by species in iteration
  df_ddP.sb <- df_ddP[df_ddP$smpldscr04==sp,] 
  # get copies for only std
  expstd <- df_ddP.sb$stdlvl[df_ddP.sb$smpldscr01=="STD"]
  obsstd <- df_ddP.sb$ccpuL[df_ddP.sb$smpldscr01=="STD"]
  # plot the two std dilutin series from ddPCR and qPCR
  plot(log10(expstd),log10(obsstd))
  # make it a data frame
  df_diff01 <- data.frame(x=obsstd,y=expstd)
  # make it a data frame
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
  # make linear model for the standard dilution series
  lmd2 <- lm(formula=y~x, data=df_diff02.1)
  # get the intercept of the linear model
  intcp2 <- lmd2$coefficients[1]
  # copy the intercept into a factor
  f2 <- 1/(10^intcp2)
  #ddlod <- min(expstd[expstd>=1E1])
  ddlod <- min(obsstd[obsstd>0])
  # collect values for each species
  lst_intcp2[[sp]] <- intcp2 
  lst_f2[[sp]] <- f2 
  lst_ddlod[[sp]] <- ddlod 
  # end iteration over subsets
}
# make the lists data frames instead of lists
df_intcp2 <- as.data.frame(do.call(rbind,lst_intcp2))
df_f2 <- as.data.frame(do.call(rbind,lst_f2))
df_ddlod <- as.data.frame(do.call(rbind,lst_ddlod))
# put all data frames into list
lst_dfs <- list(df_intcp2, df_f2, df_ddlod)
# bind data frames together
df_ifq <- do.call(cbind, lst_dfs)
# change column names
colnames(df_ifq) <- c("intcp2","f2","ddlod")
# get the species abbreviation
df_ifq$spcAbbr <- row.names(df_ifq)
#unique(df_e02$speciesabbr) %in% unique(df_ifq$spcAbbr)
# match to get modification factor 'f2' and 'ddlod' and 'intercept'
df_e02$f2 <- df_ifq$f2[match(df_e02$speciesabbr,df_ifq$spcAbbr)]
df_e02$ddlod <- df_ifq$ddlod[match(df_e02$speciesabbr,df_ifq$spcAbbr)]
df_e02$intcp2 <- df_ifq$intcp2[match(df_e02$speciesabbr,df_ifq$spcAbbr)]
# modify the copy count in the qPCR with the species corresponding modification factor
df_e02$cppuLm <- df_e02$cppuL*df_e02$f2

# get factor for adjusting qPCR LOD and qPCR LOQ
df_qP.lod.loq$f2 <- df_ifq$f2[match(df_qP.lod.loq$spcabbr,df_ifq$spcAbbr)]
#adjust the LOD and LOQ
df_qP.lod.loq$LOD_f2<- df_qP.lod.loq$f2*df_qP.lod.loq$LOD_mean
df_qP.lod.loq$LOQ_f2<- df_qP.lod.loq$f2*df_qP.lod.loq$LOQ_mean

df_e02.tmp <- df_e02[df_e02$speciesabbr=="Karmik",]
#View(df_e02.tmp)
# use dplyr to select columns and group by these columns, and # then tak the mean, the stddev, and stderr
df_e03 <- df_e02 %>%
  dplyr::select(speciesabbr, Replicate, cppuLm,smpltp,WellType ) %>% 
  dplyr::group_by(speciesabbr,Replicate,smpltp) %>%
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
# get the template volume
df_ddP$templvol <- as.numeric(gsub("voltempl_([0-9]{1})uL","\\1",df_ddP$smpldscr02))

# copy the column with standard dilution level
#df_ddP$stdlvl <- df_ddP$Sample.description.3
df_ddP$stdlvl[is.na(df_ddP$stdlvl)] <- 0

df_ddP$dsc2std <- df_ddP$stdlvl
#dsc2 <- as.character(log10(df_ddP$Sample.description.3))
#df_ddP$dsc2.2 <- gsub("^","STD5E",dsc2)
df_ddP$PoissonConfMin <- as.numeric(df_ddP$PoissonConfMin) 
df_ddP$PoissonConfMax <- as.numeric(df_ddP$PoissonConfMax) 
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
#df_e03[df_e03$qP.speciesabbr=="Karmik",]
# copy data frame
df_e04 <- df_e03
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
dfg$Sample.description.1 <- gsub("std","",dfg$Sample.description.1)
# replace the column with smpl.No for the rows in 'Sample.description.1 that match 'STD'   
df_ddP$smpl.No[df_ddP$Sample.description.1=="STD"] <- dfg$Sample.description.1

# I used 5 uL of template so the confidence intervals should be multiplied with 5
df_ddP$ccpuL <- df_ddP$Conc.copies.µL.*5
#df_ddP$ccpuL <- df_ddP$Copies.20µLWell/4
# substitute in the negative extraction control
df_ddP$smpl.No[grepl("NEK",df_ddP$smpl.No)] <- gsub("Feb","feb", df_ddP$smpl.No[grepl("NEK",df_ddP$smpl.No)])

# find sampletypes included in ddPCR setup not included in qPCR setup
# and addthem to vector
extsmplTp <- df_ddP$smpl.No[!df_ddP$smpl.No %in% df_e04$qP.smpltp]
# get the 'tidyverse' functions
library(tidyverse)
# add rows to the data frame: https://stackoverflow.com/questions/28467068/how-to-add-a-row-to-a-data-frame-in-r
df_e04.1 <- as.data.frame(df_e04) %>% add_row("qP.smpltp" = extsmplTp)
# make it a tibble again
df_e04 <- as_tibble(df_e04.1)
# modify the sample description for the std levels
stdlvl.dd <- df_ddP$smpldscr03[df_ddP$Sample.description.6=="STD"]
stdlvl.dd <- gsub("^1(E)","std\\1",stdlvl.dd)
df_ddP$Sample.description.6[df_ddP$Sample.description.6=="STD"] <-  stdlvl.dd
# paste together species abbreviation and sample name
df_e04$spcAbbr.smplNo <- paste0(df_e04$qP.speciesabbr,"_", df_e04$qP.smpltp)
df_ddP$spcAbbr.smplNo <- paste0(df_ddP$smpldscr04,"_", df_ddP$Sample.description.6)
# match between data frame to combine qPCR and ddPCR results
df_e04$dP.ccp <- df_ddP$ccpuL[match(df_e04$spcAbbr.smplNo,df_ddP$spcAbbr.smplNo)]
# replace NA in the columns that has values
df_e04 %>% replace_na(list(qP.mean=0,
                           qP.sd=0,
                           qP.se=0,
                           dP.ccp=0))

# copy the standard deviation for the qPCR results to have a minimum and a maximum
# to match the ddPCR 
df_e04$qP.sdmn <- df_e04$qP.mean-df_e04$qP.sd
df_e04$qP.sdmx <- df_e04$qP.mean+df_e04$qP.sd
# match to get confidence intervals for ddPCR evaluations
df_e04$dP.pcmn <- df_ddP$PoissonConfMin[match(df_e04$spcAbbr.smplNo,df_ddP$spcAbbr.smplNo)]
df_e04$dP.pcmx <- df_ddP$PoissonConfMax[match(df_e04$spcAbbr.smplNo,df_ddP$spcAbbr.smplNo)]
# make 'PoissonConfidenceMax68' numeric
df_ddP$PoissonConfidenceMin68 <- as.numeric(df_ddP$PoissonConfidenceMin68)
df_ddP$PoissonConfidenceMax68 <- as.numeric(df_ddP$PoissonConfidenceMax68)
# I used 5 uL of template so the confidence intervals should be multiplied with 5
df_ddP$PCmn68t5 <- df_ddP$PoissonConfidenceMin68*5
df_ddP$PCmx68t5 <- df_ddP$PoissonConfidenceMax68*5
# match back to get upper and lower confidence ingterval
df_e04$dP.pcmn <- df_ddP$PCmn68t5[match(df_e04$spcAbbr.smplNo,df_ddP$spcAbbr.smplNo)]
df_e04$dP.pcmx <- df_ddP$PCmx68t5[match(df_e04$spcAbbr.smplNo,df_ddP$spcAbbr.smplNo)]
#colnames(df_e04)
# replace NA with zeroes in multiple columns
df_e04 <- df_e04 %>% 
  mutate_at(c('dP.ccp','dP.pcmn','dP.pcmx'), ~replace_na(.,0))
# remove a column
df_e04$spcAbbr.smplNo <- NULL
# use select in dplyr to keep specified columns. dplyr requires a grouping variable
df_e04 <- df_e04 %>%  dplyr::group_by(qP.speciesabbr, qP.smpltp) %>% 
  dplyr::select(qP.speciesabbr,qP.smpltp,qP.mean,dP.ccp,qP.sdmn,qP.sdmx,dP.pcmn,dP.pcmx) 
# rename column headers
colnames(df_e04) <- c("speciesabbr","smplNm","qP.mcp",
                      "dP.mcp","qP.sdmn","qP.sdmx","dP.sdmn","dP.sdmx")
# use gather in tidyr package to  
df_e05 <- df_e04 %>% tidyr::gather(key = "cattype", value = "cpcnt", -c(speciesabbr,smplNm))
# exclude the NAs from the data frame
df_e05 <- df_e05[!is.na(df_e05$speciesabbr),]
# get unique species names
sps <- unique(df_e05$speciesabbr)
# modify in the list of the csv input files
lst_inf03 <- gsub(wd00_01,"",lst_inf02)
lst_inf03 <- gsub("/","",lst_inf03)
# substitute the wrong species abbreviation
lst_inf03 <- gsub("Myaara","Myaare",lst_inf03)
# use grepl to match vector in vector, which gives a matrix
mtx_mtch <- sapply(sps,grepl, lst_inf03)
# get only TRUE matches in this vector
# https://stackoverflow.com/questions/9505849/r-how-to-get-row-and-column-names-of-the-true-elements-of-a-matrix
spsmtch <- apply(
  mtx_mtch, 1, 
  function(u) paste( names(which(u)), collapse="," ) 
)
# combine in to a data frame
df_fls_spNm <- as.data.frame(cbind(lst_inf03,spsmtch))
# change colnames
colnames(df_fls_spNm) <- c("flNm","spcNm")
# get the unique species names in the main data frame
spcNms <- unique(df_e05$speciesabbr)
# order the species names alphabetically
spcNms <- spcNms[order(spcNms)]
# make an empty list to add dataframes to
lst_dfe <- list()
# make a count number 
i <- 1
# iterate over species names 
for (spcNm in spcNms)
  
{print(spcNm)
  #}
  df_e05.sb <- df_e05[df_e05$speciesabbr==spcNm,]
  # use separate to have two colunms, one for the machine type, and one for the category of
  # counts
  df_e05.sb <- tidyr::separate(data = df_e05.sb, 
                               col = cattype, 
                               into = c("machine", "categ"), sep = "\\.")
  # get limit of detection
  ddlod <- df_ddlod$V1[match(spcNm,rownames(df_ddlod))]
  # get the 'std' samples, among the 'qPCR' machine data, 
  # among the mean copy count 'mcp', and get the lowest
  # above zero
  df_e05.sb.std.mcp  <- df_e05.sb[(grepl("std",df_e05.sb$smplNm) & 
                                     grepl("q",df_e05.sb$machine) & 
                                     grepl("mcp",df_e05.sb$categ)),]
  # get the qlod
  qlod <- min(df_e05.sb.std.mcp$cpcnt[df_e05.sb.std.mcp$cpcnt>0])
  qloq.unmf <- df_qPCRm04$smpltp.dil[df_qPCRm04$speciesabbr==spcNm]
  sp.f2 <- df_ifq$f2[df_ifq$spcAbbr==spcNm]
  df_e05.sb.std.mcp$org.dil.lvl <- as.numeric(gsub("std","3",df_e05.sb.std.mcp$smplNm))
  qloq <- df_e05.sb.std.mcp$cpcnt[(df_e05.sb.std.mcp$org.dil.lvl==qloq.unmf)]
  # in case there is no limit then use zero instead
  if (is.na(ddlod)){
    ddlod <- 0
  }
  
  if (purrr::is_empty(qloq)){
    qloq <- 0
  }
  
  # match to get csv filename that matches the species abbreviation
  flNm <- df_fls_spNm$flNm[spcNm==df_fls_spNm$spcNm]
  # spread the data frame again based on the sub category under the machine
  df_e05.sb <- df_e05.sb %>% tidyr::spread(key = "categ", value = "cpcnt")
  # substitute to replace abbreviated machine names with long machine names
  df_e05.sb$machine <- gsub("dP","ddPCR",df_e05.sb$machine)
  df_e05.sb$machine <- gsub("qP","qPCR",df_e05.sb$machine)
  # make a list with colors for the sampling categories
  clfH <- rep(c("grey54","white"),length(unique(df_e05.sb$smplNm)))
  # count the number of elements
  ncl <- length(clfH)
  clfH[(ncl-63):(ncl)] <- rep(c("orchid4","orchid1"),63)
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
  listSmplNm <- unique(df_e05.sb$smplNm)
  # exclude control names from the list of unique sample ids
  listSmplNm <-    listSmplNm[!listSmplNm %in% controlnames]
  # add control names back to the list of unique sample ids (at the end), now in the order we want
  listSmplNm <-    c(listSmplNm, controlnames) 
  # convert smplNm to a factor with the levels (order) given by the list we have made
  df_e05.sb$smplNm <- factor(df_e05.sb$smplNm, levels= listSmplNm)
  # remember that for ggplot() we should now use “fill=smplNm” instead of “fill=as.factor(smplNm)”
  # https://stackoverflow.com/questions/41047939/ggplot-barplot-how-to-display-small-positive-numbers-with-log-scaled-y-axis
  pd = position=position_dodge(.5)
  plt07 <- ggplot(data=df_e05.sb, aes(x=smplNm,y=mcp)) +
    theme_bw() +
    theme_classic() +
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
    geom_errorbar(aes(ymin=sdmn, ymax=sdmx, colour=machine), 
                  width=0.9, lty=1,position=pd, size=0.8) +
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
    geom_hline(yintercept=ddlod, color="springgreen4",lty=2) +
    #add a text label to the horizontal line
    geom_text(aes(x = (as.factor(smplNm)[90]), y = ddlod-(ddlod/5)), size=3.4,
              label = "LOQ-STD-DIL-SER", color = "springgreen4") +
    geom_hline(yintercept=qloq, color="#E69F00",lty=2) +
    geom_text(aes(x = (as.factor(smplNm)[90]), y = qlod-(qlod/5)), size=3.4,
              label = "qLOD", color = "#E69F00") +
    geom_hline(yintercept=qlod, color="#E69F00",lty=1) +
    geom_text(aes(x = (as.factor(smplNm)[90]), y = qloq-(qloq/5)), size=3.4,
              label = "qLOQ", color = "#E69F00") +
    theme(axis.text.x = element_text(angle = 90, hjust=1))   +
    ylab("eDNA molecules/ uL extraktion from filtersample") + 
    xlab("sample number") +
    annotation_logticks(sides="l") +
    #coord_flip()+
    scale_color_manual(values=c("springgreen4", "#E69F00"))
  # change the plot header
  #plt07  <- plt07 + labs(title = paste0("eDNA detekteret med ",spcNm,"-assay i 2021 prøver"))#,
  plt07  <- plt07 + labs(title = paste0("eDNA detected using ",spcNm,"-assay in 2021"))#,
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  plt07 <- plt07 + labs(color='machine')
  plt07 <- plt07 + labs(fill='machine')
  plt07 <- plt07 + labs(shape='machine')
  
  #plt07 
  # modify the input file name
  inf02_1 <- gsub("\\.csv","",flNm)
  # get the count number
  # and pad with zeros to two characters
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pltnmb <-stringr::str_pad(i, 2, pad = "0")
  # add one increment to the count number
  i <- i+1
  #make filename to save plot to save the plot to
  fgNm07 <- paste0(wd00_01,"/Fig02_v02_plot",pltnmb,"_",inf02_1,".png")
  # save the plot as png file if 'bSaveFigures' is set to be TRUE
  bSaveFigures <- TRUE
  if(bSaveFigures==T){
    ggsave(plt07,file=fgNm07,#width=210,height=297,
           width=297,height=210,
           units="mm",dpi=600)
  }
  
  #end iteration over species names
  #}
  # colnames(df_e05)
  #_______________________________________________________________________________
  # example with tidyr
  #_______________________________________________________________________________
  library(tidyr)
  # #_______________________________________________________________________________
  # use tidyr example above
  df_e06.sb  <- df_e05.sb %>%
    tidyr::pivot_wider(names_from="machine",
                       values_from=c("mcp","sdmn","sdmx"))
  # add the data frame to the list of dataframes
  lst_dfe[[i]] <- df_e06.sb
  # get the library for making plots
  library(ggplot2)
  # add a column for the sample type
  df_e06.sb$smplTp <- "MST"
  df_e06.sb$smplTp[grepl("std",df_e06.sb$smplNm)] <- "STD"
  plt08 <- ggplot(data=df_e06.sb, aes(x=mcp_qPCR ,y=mcp_ddPCR) ) +
    theme_classic() +
    geom_point(aes(fill=smplTp, shape=smplTp), size=3) + 
    scale_shape_manual(values = c(21, 24)) +
    scale_fill_manual(values=c("red", "pink")) +
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
    #xlab("antal eDNA molekyler /uL i qPCR") + ylab("antal eDNA molekyler /uL i ddPCR") +
    xlab("eDNA molecules/uL in qPCR") + ylab("eDNA molecules/uL in ddPCR") +
    annotation_logticks(sides="bl") +
    geom_abline(intercept = 0, slope = 1, color="blue") +
    # add line for LOQ and/ or LOD
    #geom_hline(yintercept=ddlod, color="springgreen4") +
    geom_segment(aes(x=0,xend=ddlod,y=ddlod,yend=ddlod),color="springgreen4",lty=2) +
    geom_text(aes(x = 1E-1, y = ddlod+(ddlod/2)), size=3.4,
              label = "LOQ-STD-DIL-SER", color = "springgreen4") +
    #geom_vline(xintercept=qloq, color="#E69F00",lty=2) +
    geom_segment(aes(x=qlod,xend=qlod,y=0,yend=qlod),color="#E69F00",lty=1) +
    geom_text(aes(y = 0.5E-1, x = qlod+(qlod/2)), size=3.4,
              angle=90, label = "qLOD", color = "#E69F00") +
    geom_segment(aes(x=qloq,xend=qloq,y=0,yend=qloq),color="#E69F00",lty=2) +
    #geom_vline(xintercept=qlod, color="#E69F00",lty=1) +
    geom_text(aes(y = 1E-1, x = qloq+(qloq/2)), size=3.4,
              angle=90,label = "qLOQ", color = "#E69F00") +
    labs(color='type') +
    labs(fill='type') +
    labs(shape='type') +
    labs(title = paste0("eDNA levels for ",spcNm,"-assay fra 2021"))#,
  
  #plt08
  fgNm08 <- paste0(wd00_01,"/Fig03_v02_plot",pltnmb,"_",inf02_1,".png")
  # save the plot as png file if 'bSaveFigures' is set to be TRUE
  bSaveFigures <- TRUE
  if(bSaveFigures==T){
    ggsave(plt08,file=fgNm08,#width=210,height=297,
           width=297,height=210,
           units="mm",dpi=600)
  }
  
  
  # #_______________________________________________________________________________
  # mchty <- c("ddPCR","qPCR") 
  # MSTsmpl071 <- df_e05$mcp[grepl("MST2021-071",df_e05$smplNm)]
  # MSTsmpl058 <- df_e05$mcp[grepl("MST2021-058",df_e05$smplNm)]
  # MSTsmpl085 <- df_e05$mcp[grepl("MST2021-085",df_e05$smplNm)]
  # df_MSTres <- as.data.frame(rbind(mchty,MSTsmpl071,MSTsmpl058,MSTsmpl085))
  # colnames(df_MSTres) <- df_MSTres[1,]
  # df_MSTres <- df_MSTres[-1,]
  
  # ------------- plot Combined figure -------------
  # USe patchwork to make combined figure
  library(patchwork)
  # set a variable to TRUE to determine whether to save figures
  bSaveFigures <- T
  #getwd()
  #define a filename to save to
  fnm02 <- paste0(wd00_01,"/Fig04_v02_plot",pltnmb,"_",inf02_1,".png")
  #see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
  # on how to arrange plots in patchwork
  
  plt07t <- plt07 + labs(title = "A")#,
  plt08t <- plt08 + labs(title = "B")#,
  p <-  plt07t  +
    plt08t +
    
    plot_layout(nrow=1,byrow=T) + #xlab(xlabel) +
    plot_layout(guides = "collect") #+
  #patchwork::plot_annotation(caption="Fig0023_v01_Psefar") #& theme(legend.position = "bottom")
  #p
  if(bSaveFigures==T){
    ggsave(p,file=fnm02,width=297*1.2,height=210*0.8,
           units="mm",dpi=300)
  }
  
  # ------------- plot Combined figure -------------
  # USe patchwork to make combined figure
  library(patchwork)
  # set a variable to TRUE to determine whether to save figures
  bSaveFigures <- T
  #getwd()
  #define a filename to save to
  fnm02 <- paste0(wd00_01,"/Fig05_v02_plot",pltnmb,"_",inf02_1,".png")
  #see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
  # on how to arrange plots in patchwork
  
  plt07t <- plt07 + labs(title = "A")#,
  plt08t <- plt08 + labs(title = "B")#,
  p <-  plt07t  +
    plt08t +
    
    plot_layout(nrow=1,
                ncol=2,byrow=T) + #xlab(xlabel) +
    plot_layout(guides = "collect") #+
  #patchwork::plot_annotation(caption="Fig0023_v01_Psefar") #& theme(legend.position = "bottom")
  #p
  if(bSaveFigures==T){
    ggsave(p,file=fnm02,
           #width=297*0.8,height=210*1.2,
           #height=297,width=210,
           width=297,height=210*0.8,
           units="mm",dpi=300)
  }
  
  #end iteration over species names
}

#_______________________________________________________________________________

# make the lists data frames instead of lists
dfe.1 <- as.data.frame(do.call(rbind,lst_dfe))
df_e07 <- dfe.1
# save the results in data frame - that can be used later on for analysis
# get limit of detection
df_e07$ddlod <- df_ddlod$V1[match(df_e07$speciesabbr,rownames(df_ddlod))]
# above zero
df_e06.1  <- df_e05[(grepl("std",df_e05$smplNm) & 
                       grepl("q",df_e05$cattype)),]
# get the qlod
df_e06.1.abz <- df_e06.1[df_e06.1$cpcnt>0,]
df_e06.1.qlod <- df_e06.1.abz %>% 
  dplyr::select(speciesabbr, cpcnt ,smplNm) %>% 
  dplyr::group_by(speciesabbr) %>%
  dplyr::summarise_each(min, c(cpcnt))
# match to get qlod
df_e07$qlod <- df_e06.1.qlod$cpcnt[match(df_e07$speciesabbr,df_e06.1.qlod$speciesabbr)]
# match to get qloq
df_e06.1$qloq.unmf <- 
  df_qPCRm04$smpltp.dil[match(df_e06.1$speciesabbr,df_qPCRm04$speciesabbr)]
# get original dilution level
df_e06.1$org.dil.lvl <- as.numeric(gsub("std","3",df_e06.1$smplNm))
df_e05.std.qP.LOQ <- df_e06.1[df_e06.1$org.dil.lvl==df_e06.1$qloq.unmf,]
df_e05.std.qP.LOQ <- df_e05.std.qP.LOQ[df_e05.std.qP.LOQ$cattype=="qP.mcp",]
# get the modified copy count for qPCR LOQ 
df_e07$qloq <- df_e05.std.qP.LOQ$cpcnt[match(df_e07$speciesabbr,df_e05.std.qP.LOQ$speciesabbr)]
# get modification factor
df_e07$sp.f2 <- df_ifq$f2[match(df_e07$speciesabbr,df_ifq$spcAbbr)]

df_e08 <- df_e07
# make a filename to save the resulting
csvflNm <- paste0("results_from_all_assays_in_ddPCR_v02.csv")
# save the resulting dat frame as a csv file
csvflNm <- paste0(wd00_01,"/",csvflNm)

write_csv(df_e08,file=csvflNm)
