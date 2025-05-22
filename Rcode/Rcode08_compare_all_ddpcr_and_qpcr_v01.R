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
inpf01  <- "outfile02_merged_qpcr_csvfls_from_MONIS_MST_2021_MST_samples.txt"
#paste path and input file together
if01wp <- paste(wd02.1,"/",inpf01,sep="")
#read in csv with qPCR data from MxPro
df_e01 <- read.csv(if01wp,header=T,sep=";")

# copy data frame
df_e02 <-  df_e01
# make the column numeric
df_e02$Quantitycopies <- as.numeric(df_e02$Quantitycopies)
# read in an xlsx file with all primer assays listed 
# for each species anmd the abbreviations
inf04 <- "list_of_specific_assays_MONIS6.xlsx"
wd00.1_inf04 <- paste0(wd00.1,"/",inf04 )
library(xlsx)
df_dtc_asss <- xlsx::read.xlsx(wd00.1_inf04,1)

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
# View(df_e02.1)
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
# then the column 'Conc.copies.µL.' should 
df_ddP$ddcpuL <-df_ddP$Conc.copies.µL.
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
# modify the copy count in the qPCR with the species corresponding 
# modification factor
df_e02$cppuLm <- df_e02$cppuL*df_e02$f2

# get factor for adjusting qPCR LOD and qPCR LOQ
df_qP.lod.loq$f2 <- df_ifq$f2[match(df_qP.lod.loq$spcabbr,df_ifq$spcAbbr)]
#adjust the LOD and LOQ
df_qP.lod.loq$LOD_f2<- df_qP.lod.loq$f2*df_qP.lod.loq$LOD_mean
df_qP.lod.loq$LOQ_f2<- df_qP.lod.loq$f2*df_qP.lod.loq$LOQ_mean

df_e02.tmp <- df_e02[df_e02$speciesabbr=="Karmik",]
# make a column that holds the entered quantities for the std dilution
# series, and the quantities inferred for the wells.
df_e02$cppuLq <- df_e02$Quantitycopies
# use dplyr to select columns and group by these columns, and 
# then take the mean, the stddev, and stderr
# https://stackoverflow.com/questions/77217041/summarize-across-with-a-specific-formula-and-missing-values

df_e03 <- df_e02 %>%
  dplyr::select(speciesabbr, Replicate, cppuLq,cppuL,cppuLm,smpltp,WellType ) %>% 
  dplyr::group_by(speciesabbr,Replicate,smpltp) %>%
  summarise(across(starts_with("cppuL"), list(mean = mean,
# for some reason the 'se' would not work without the '~' in front.                                              
                                  se =  ~ sd(.x, na.rm = TRUE)/sqrt(length(.x)),
                                              sd = sd),
                                              .names = "{.col}.{.fn}"))
  #dplyr::summarise(dplyr::across(c(cppuLm,cppuL), list(mean=mean, sd=sd, se=sd/sqrt(n()))))
  #dplyr::summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), c(cppuLm,cppuL))

# First, create a list of all column names and set to 0
lst_df_e03 <- setNames(lapply(vector("list", ncol(df_e03)),
                          function(x) x <- 0), names(df_e03))

# Now use that list in tidyr::replace_na 
df_e03 <- df_e03 %>% replace_na(lst_df_e03)
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

# get the 'tidyverse' functions
library(tidyverse)
# make it a tibble again
df_e04 <- as_tibble(df_e04)
#View(df_e04)
# modify the sample description for the std levels
stdlvl.dd <- df_ddP$smpldscr03[df_ddP$Sample.description.6=="STD"]
stdlvl.dd <- gsub("^1(E)","std\\1",stdlvl.dd)
df_ddP$Sample.description.6[df_ddP$Sample.description.6=="STD"] <-  stdlvl.dd
# paste together species abbreviation and sample name
df_e04$spcAbbr.smplNo <- paste0(df_e04$qP.speciesabbr,"_", df_e04$qP.smpltp)
df_ddP$spcAbbr.smplNo <- paste0(df_ddP$smpldscr04,"_", df_ddP$Sample.description.6)

# use the 'Sample.description.6' as the smpl number colunm
df_ddP$smpl.No <- df_ddP$Sample.description.6
# find sampletypes included in ddPCR setup not included in qPCR setup
# and addthem to vector
extsmplTp <- df_ddP$smpl.No[!df_ddP$smpl.No %in% df_e04$qP.smpltp]

# add rows to the data frame: https://stackoverflow.com/questions/28467068/how-to-add-a-row-to-a-data-frame-in-r
df_e04 <- as.data.frame(df_e04) %>% add_row("qP.smpltp" = extsmplTp)


# match between data frame to combine qPCR and ddPCR results
df_e04$dP.ccp <- df_ddP$ccpuL[match(df_e04$spcAbbr.smplNo,df_ddP$spcAbbr.smplNo)]

# copy the standard deviation for the qPCR results to have a minimum and a maximum
# to match the ddPCR 
df_e04$qP.cppuL.sdmn <- df_e04$qP.cppuL.mean-df_e04$qP.cppuL.sd
df_e04$qP.cppuL.sdmx <- df_e04$qP.cppuL.mean+df_e04$qP.cppuL.sd

df_e04$qP.cppuLq.sdmn <- df_e04$qP.cppuLq.mean-df_e04$qP.cppuLq.sd
df_e04$qP.cppuLq.sdmx <- df_e04$qP.cppuLq.mean+df_e04$qP.cppuLq.sd

df_e04$qP.cppuLm.sdmn <- df_e04$qP.cppuLm.mean-df_e04$qP.cppuLm.sd
df_e04$qP.cppuLm.sdmx <- df_e04$qP.cppuLm.mean+df_e04$qP.cppuLm.sd


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
  dplyr::select(qP.speciesabbr,
                qP.smpltp,
                qP.cppuL.mean,
                qP.cppuLm.mean,
                qP.cppuLq.mean,
                dP.ccp,
                qP.cppuL.sdmn,
                qP.cppuLm.sdmn,
                qP.cppuLq.sdmn,
                qP.cppuL.sdmx,
                qP.cppuLm.sdmx,
                qP.cppuLq.sdmx,
                dP.pcmn,
                dP.pcmx) 
# # rename column headers
colnames(df_e04)[grepl("speciesabbr",colnames(df_e04))] <- c("speciesabbr")
colnames(df_e04)[grepl("smpl",colnames(df_e04))] <- c("smplNm")
# colnames(df_e04) <- c("speciesabbr","smplNm","qP.mcp",
#                       "dP.mcp","qP.sdmn","qP.sdmx","dP.sdmn","dP.sdmx")
# use gather in tidyr package to  
df_e05 <- df_e04 %>% tidyr::gather(key = "cattype", 
                                   value = "cpcnt", -c(speciesabbr,smplNm))
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



# make a function, that can be used in the dplyr pipe
div5 <- function(x, na.rm = FALSE) (x/5)
# load the thedplyt library
library(dplyr)
# for all columns that starts with "dP." apply the function
# that divides by 5. 
# see this example : https://forum.posit.co/t/applying-mutate-across-starts-with-parse-number-all-together/126758/4
df_e04a <- df_e04 %>%
  dplyr::mutate_at(dplyr::vars(dplyr::starts_with("dP.")), div5)
# then group by the 'speciesabbr', 'smplNm' columns
# and the only select specific columns
df_e04a <- df_e04a %>%  dplyr::group_by(speciesabbr,
                                        smplNm) %>% 
  dplyr::select(speciesabbr,
                smplNm,
                qP.cppuL.mean,
                #qP.cppuLm.mean,
                #qP.cppuLq.mean,
                dP.ccp,
                qP.cppuL.sdmn,
                #qP.cppuLm.sdmn,
                #qP.cppuLq.sdmn,
                qP.cppuL.sdmx,
                #qP.cppuLm.sdmx,
                #qP.cppuLq.sdmx,
                dP.pcmn,
                dP.pcmx) 

# change the column headers
clNmse04a <- colnames(df_e04a)
clNmse04a[grepl("dP.ccp",clNmse04a)] <- c("dP.cppuL.mean")
clNmse04a[grepl("dP.pc",clNmse04a)] <- c("dP.cppuL.sdmn", "dP.cppuL.sdmx")
clNmse04a[grepl("mean",clNmse04a)] <- gsub("mean","mlcnt",clNmse04a[grepl("mean",clNmse04a)])
# use 'sub' instead og 'gsub' to replace only first occurence
clNmse04a <- sub("\\.","_",clNmse04a)
# use the new column names
colnames(df_e04a) <- clNmse04a
# I did not get the 'Neomel', the 'Rhihar' and the 'Oncgor' 
# assays working in this setup, and because of this I am excluding them here
df_e04a <- df_e04a[!(df_e04a$speciesabbr=="Neomel"),]
df_e04a <- df_e04a[!(df_e04a$speciesabbr=="Rhihar"),]
df_e04a <- df_e04a[!(df_e04a$speciesabbr=="Oncgor"),]

# copy the column with species abbreviation names in to 
# a new column
df_qP.lod.loq$speciesabbr <- df_qP.lod.loq$spcabbr
# match to get ddPCR lod
df_qP.lod.loq$ddlod <- df_ddlod$V1[match(df_qP.lod.loq$speciesabbr,rownames(df_ddlod))]
# use dplyr::left_join to match between data frames 
# see this example : https://stackoverflow.com/questions/32209396/how-to-select-all-columns-in-dplyr-sql
df_e04a <- dplyr::left_join(df_e04a,
                         # select all columns to include           
                         df_qP.lod.loq %>% dplyr::select(everything()),
                         by = "speciesabbr")

# substitute to change in the column names
clNme04a <- colnames(df_e04a)
clNme04a <- gsub("_","",clNme04a)
clNme04a[grepl("LO",clNme04a)] <- paste0("qPcppuL.",clNme04a[grepl("LO",clNme04a)])
clNme04a[grepl("ddlod",clNme04a)] <- paste0("dPcppuL.",clNme04a[grepl("ddlod",clNme04a)])
colnames(df_e04a) <- clNme04a

# load package libraries
library(tibble)
library(tidyr)
library(dplyr)
colnames(df_e04a)

colnames(df_e04a)[grepl("dP",colnames(df_e04a))]

# use tidyr::pivot_longer
df_e04b <- df_e04a %>%
  # do not select unneeded columns
  dplyr::select(-c(spcabbr,f2)) %>%
  tidyr::pivot_longer(.,
                      -c(speciesabbr,
                         smplNm), 
                      #cols = starts_with(c("qP_cppuL","dP_cppuL")),
                      #names_sep = "\\.",
                      #names_sep = "_",
                      names_to = c('mch','vll'), 
                      values_to = "molcnt",
                      #names_pattern = '(.*)\\.(t\\d+)') %>%
                      names_pattern = '(.*)\\.(.*)') %>%
  dplyr::group_by(smplNm)
df_e04b$molcnt <- as.integer(df_e04b$molcnt)
# reorder the data frame
df_e04b <- df_e04b %>% dplyr::arrange(speciesabbr,smplNm, mch,vll)

# rearrange the longer version in to a wider version
# se:  https://stackoverflow.com/questions/61367186/pivot-longer-into-multiple-columns
df_e04c <- df_e04b %>% pivot_wider(., id_cols = c(speciesabbr,
                                     smplNm, mch), 
                         names_from = vll, 
                         values_from = molcnt, 
                         names_repair = "check_unique")
# change in the column with the machine
df_e04c$mch <- gsub("dPcppuL","ddPCR",df_e04c$mch)
df_e04c$mch <- gsub("qPcppuL","qPCR",df_e04c$mch)
# use the new species name assigned to 'Magallana	gigas'
df_e04c$speciesabbr <- gsub("Cragig","Maggig",df_e04c$speciesabbr)
# use the data frame with the detection assays to get the full
# length Latin Genus and Species names
df_e04c$GnNm <-  df_dtc_asss$Genus[match(df_e04c$speciesabbr,df_dtc_asss$AbbrvNm)]
df_e04c$SpNm <-  df_dtc_asss$Species[match(df_e04c$speciesabbr,df_dtc_asss$AbbrvNm)]
df_e04c$GnNmLt <- substr(df_e04c$GnNm, 1, 1)
df_e04c$speciesNm <- paste0(df_e04c$GnNmLt,". ",df_e04c$SpNm)
# make a list with colors for the sampling categories
clfH <- rep(c("grey54","white"),length(unique(df_e04c$smplNm)))
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
# remember that for ggplot() we should now use “fill=smplNm” instead of “fill=as.factor(smplNm)”
# https://stackoverflow.com/questions/41047939/ggplot-barplot-how-to-display-small-positive-numbers-with-log-scaled-y-axis
pd = position=position_dodge(.5)
df_e04d <- df_e04c[!is.na(df_e04c$LODmean),]

#df_e05
plt07 <- ggplot(data=df_e04c, aes(x=smplNm,y=mlcnt )) +
  theme_bw() +
  theme_classic() +
  # justify to the left margin
  theme_minimal() +
  #https://stackoverflow.com/questions/62009919/facet-title-alignment-using-facet-wrap-in-ggplot2
  theme(strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
        panel.background = element_rect(fill = "grey99",
                                        color = "white")) +
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
  #theme(panel.background = element_rect(fill = "transparent")) +
  # https://stackoverflow.com/questions/2678141/how-can-i-suppress-the-vertical-gridlines-in-a-ggplot2-plot
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() #,
    # explicitly set the horizontal lines (or they will disappear too)
    #panel.grid.major.y = element_line( size=.1, color="black" ) 
  ) +
  theme(panel.background = element_rect(fill = "white")) +
  # add points and error bars
  geom_errorbar(aes(ymin=sdmn, ymax=sdmx, colour=mch), 
                width=0.9, lty=1,position=pd, linewidth=0.8) +
  geom_point(aes(colour=mch), position=pd, size=2) +
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
  
  # geom_hline(yintercept=df_e04c$LODmean, color="springgreen4",lty=2) +
  # #add a text label to the horizontal line
  # geom_text(aes(x = (as.factor(smplNm)[90]), y = ddlod-(ddlod/5)), size=3.4,
  #           label = "LOQ-STD-DIL-SER", color = "springgreen4") +
   # geom_hline(yintercept=df_e04d$LODmean, color="#E69F00",lty=2) +
  # geom_text(aes(x = (as.factor(smplNm)[90]), y = LODmean-(LODmean/5)), size=3.4,
  #           label = "qLOD", color = "#E69F00") +
  #  geom_hline(yintercept=LODmean, color="#E69F00",lty=1) +
  # geom_text(aes(x = (as.factor(smplNm)[90]), y = qloq-(qloq/5)), size=3.4,
  #           label = "qLOQ", color = "#E69F00") +

  ylab("eDNA molecules/ uL") + 
  xlab("sample number") +
  annotation_logticks(sides="l") +
  facet_wrap(vars(speciesNm ), ncol = 2) +
  
  theme(axis.text.x = element_text(angle = 90, hjust=1))   +
  
  
  scale_color_manual(values=c("springgreen4", "#E69F00"))
#getting separate legends
plt07 <- plt07 + labs(color='machine')
plt07 <- plt07 + labs(fill='machine')
plt07 <- plt07 + labs(shape='machine')
#plt07 

#make filename to save plot to save the plot to
fgNm07 <- paste0(wd00_01,"/Fig11_scatterplot_in_colmns_v03.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- TRUE
if(bSaveFigures==T){
  ggsave(plt07,file=fgNm07,
         width=210,height=297*1.4,
         #width=297,height=210,
         #width=297,height=210*1.4,
         units="mm",dpi=300)
}
#_______________________________________________________________________________
# substitute in the category type
df_e05$cattype <- gsub("(.P)\\.(.*)\\.","\\1.\\2_",df_e05$cattype)
# make an empty list to add dataframes to
lst_dfe <- list()
# make a count number 
i <- 1
# make an empty list
lst_dfspc.e05 <- list()
# make seq of numbers for species
nSpc <- seq(1,length(spcNms),1)
# iterate over species names 
for (n in nSpc)
  {
  spcNm <- spcNms[n]
  print(spcNm)
  df_e05.sb <- df_e05[df_e05$speciesabbr==spcNm,]
  unique(df_e05.sb$cattype)
  #}
  # use separate to have two colunms, one for the machine type, and one for the category of
  # counts
  df_e05.sb <- tidyr::separate(data = df_e05.sb, 
                               col = cattype, 
                               into = c("machine", "categ"), sep = "\\.")
  # remove unneeded rows
  df_e05.sb <- df_e05.sb[!grepl("cppuLm_",df_e05.sb$categ),]
   df_e05.sb <- df_e05.sb[!grepl("cppuL_",df_e05.sb$categ),]
  # df_e05.sb <- df_e05.sb[!grepl("cppuLq_",df_e05.sb$categ),]
  # substitute
  df_e05.sb$categ <- gsub("cppuL._","cppuL_",df_e05.sb$categ)
  df_e05.sb$categ <- gsub(".*mean","mcp",df_e05.sb$categ)
  df_e05.sb$categ <- gsub("ccp","mcp",df_e05.sb$categ)
  df_e05.sb$categ <- gsub(".*mx","sdmx",df_e05.sb$categ)
  df_e05.sb$categ <- gsub(".*mn","sdmn",df_e05.sb$categ)
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
  df_e05.sb <- 
    df_e05.sb %>% tidyr::spread(key = "categ", value = "cpcnt")
  # Spread and gather are complements
  # substitute to replace abbreviated machine names with long machine names
  df_e05.sb$machine <- gsub("dP","ddPCR",df_e05.sb$machine)
  df_e05.sb$machine <- gsub("qP","qPCR",df_e05.sb$machine)
  # get unique sample ids from the dataframe used for plotting
  listSmplNm <- unique(df_e05.sb$smplNm)
  # exclude control names from the list of unique sample ids
  listSmplNm <-    listSmplNm[!listSmplNm %in% controlnames]
  # add control names back to the list of unique sample ids (at the end), now in the order we want
  listSmplNm <-    c(listSmplNm, controlnames) 
  # convert smplNm to a factor with the levels (order) given by the list we have made
  df_e05.sb$smplNm <- factor(df_e05.sb$smplNm, levels= listSmplNm)
  # collect 
  df_e05.sb$ddlod <- ddlod
  df_e05.sb$qlod <- qlod
  df_e05.sb$qloq <- qloq
  df_e05.sb$qloq.unmf <- qloq.unmf
  lst_dfspc.e05[[n]] <- df_e05.sb 
}  

#bind the rows in each list in to one data frame
df_spc.e05 <- data.table::rbindlist(lst_dfspc.e05, fill=T)
df_spc.e05 <- as.data.frame(df_spc.e05)
#View(df_spc.e05)
# only keep selected columns to make data frame that holds
# the limit of detection and quatnification
df_lo.dPqP <- df_spc.e05 %>% dplyr::select(speciesabbr,
                              ddlod,qlod,qloq,qloq.unmf)
# and only keep the unique rows
df_lo.dPqP <- df_lo.dPqP  %>% dplyr::distinct()


#_______________________________________________________________________________
  # example with tidyr
  #_______________________________________________________________________________
  library(tidyr)
  # #_______________________________________________________________________________
  # use tidyr example above
  df_e06.sb  <- df_spc.e05 %>%
    tidyr::pivot_wider(names_from="machine",
                       values_from=c("mcp","sdmn","sdmx"))
  
  # use the new species name assigned to 'Magallana	gigas'
  df_e06.sb$speciesabbr <- gsub("Cragig","Maggig",df_e06.sb$speciesabbr)
  # I did not get the 'Neomel', the 'Rhihar' and the 'Oncgor' 
  # assays working in this setup, and because of this I am excluding them here
  df_e06.sb <- df_e06.sb[!(df_e06.sb$speciesabbr=="Neomel"),]
  df_e06.sb <- df_e06.sb[!(df_e06.sb$speciesabbr=="Rhihar"),]
  df_e06.sb <- df_e06.sb[!(df_e06.sb$speciesabbr=="Oncgor"),]
  # use the data frame with the detection assays to get the full
  # length Latin Genus and Species names
  df_e06.sb$GnNm <-  df_dtc_asss$Genus[match(df_e06.sb$speciesabbr,df_dtc_asss$AbbrvNm)]
  df_e06.sb$SpNm <-  df_dtc_asss$Species[match(df_e06.sb$speciesabbr,df_dtc_asss$AbbrvNm)]
  df_e06.sb$GnNmLt <- substr(df_e06.sb$GnNm, 1, 1)
  df_e06.sb$speciesNm <- paste0(df_e06.sb$GnNmLt,". ",df_e06.sb$SpNm)
  
  # add the data frame to the list of dataframes
  #lst_dfe[[i]] <- df_e06.sb
  # get the library for making plots
  library(ggplot2)
  # add a column for the sample type
  df_e06.sb$smplTp <- "MST"
  df_e06.sb$smplTp[grepl("std",df_e06.sb$smplNm)] <- "STD"
  # # join data frames to get the lod and loq for qPCR and for ddPCR
  df_e06.sb <- df_e06.sb %>% dplyr::left_join(df_lo.dPqP,
                                  by = "speciesabbr",
                                 keep=F)
  # only keep one set of columns
  df_e06.sb <- df_e06.sb[!grepl("\\.y",colnames(df_e06.sb))]
  # substitute to remove the added suffix
  colnames(df_e06.sb) <- gsub("(.*)\\..*","\\1",colnames(df_e06.sb))
  #unique(df_e06.sb$smplTp)
  
  #View(df_e06.sb)
  plt08 <- ggplot(data=df_e06.sb, aes(x=mcp_qPCR ,y=mcp_ddPCR) ) +
    theme_classic() +
    # justify to the left margin
    theme_minimal() +
    #https://stackoverflow.com/questions/62009919/facet-title-alignment-using-facet-wrap-in-ggplot2
    theme(strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
          panel.background = element_rect(fill = "grey99",
                                          color = "white")) +
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
    # plot the points above all the other vlines
    geom_point(aes(fill=smplTp, shape=smplTp), size=2.4) + 
    scale_shape_manual(values = c(21, 24)) +
    scale_fill_manual(values=c("red", "pink")) +
    geom_errorbar(aes(ymin =  sdmn_ddPCR,ymax =  sdmx_ddPCR)) + 
    geom_errorbarh(aes(xmin =sdmn_qPCR,xmax = sdmx_qPCR)) +
    
    geom_abline(intercept = 0, slope = 1, color="blue",size=0.22,lty=2) +
    # add line for LOQ and/ or LOD
    # geom_hline(yintercept=ddlod, color="springgreen4") +
     geom_segment(aes(x=0,xend=ddlod,y=ddlod,yend=ddlod),
      color="springgreen4",lty=1,size=0.32) +
    # geom_text(aes(x = 1E-1, y = ddlod+(ddlod/2)), size=2.4,
    #           label = "LOQ-STD-DIL-SER", color = "springgreen4") +
    geom_segment( aes(x=qlod,
                     xend=qlod,y=0,
                     yend=qlod),color="#E69F00",lty=1,size=0.32) +
    # geom_text(aes(y = 0.5E-1, x = qlod+(qlod/2)), size=2.4,
    #           angle=90, label = "qLOD", color = "#E69F00") +
    geom_segment(aes(x=qloq.unmf,
                     xend=qloq.unmf,y=0,
                     yend=qloq.unmf),color="#E69F00",lty=2,size=0.32) +
    # geom_text(aes(y = 1E-1, x = qloq.unmf+(qloq.unmf/2)), size=2.4,
    #           angle=90,label = "qLOQ", color = "#E69F00") +
    # add a linear model : https://stackoverflow.com/questions/15633714/adding-a-regression-line-on-a-ggplot
    # and see: https://stackoverflow.com/questions/9613578/change-standard-error-color-for-geom-smooth
    # geom_smooth(method='lm',span = 0.9, 
    #             lty=1,size=0.3,
    #             aes(fill = smplTp,
    #                 color = smplTp)) +
    # 
    geom_smooth(method='lm',span = 0.9, 
                lty=1,size=0.3, fill="tomato3", color="red") +
  
    facet_wrap(vars(speciesNm ), ncol = 2) +
    # # https://stackoverflow.com/questions/17073772/ggplot2-legend-on-top-and-margin
    # theme(
    #   legend.margin=unit(-0.6,"cm"),
    #   #plot.background=element_rect(fill="red"),
    #   legend.position="top") +
    # guides(fill=guide_legend(title.position="top")) +
    #labs(title = paste0("eDNA levels for ",spcNm,"-assay fra 2021")) +
    labs(color='type') +
    labs(fill='type') +
    labs(shape='type') 
  #plt08
  # change background of plot
  # see : https://stackoverflow.com/questions/6736378/how-do-i-change-the-background-color-of-a-plot-made-with-ggplot2
  plt08 <- plt08 + theme(panel.background = element_rect(fill='white'))  
  
  fgNm08 <- paste0(wd00_01,"/Fig12_v01_scatterplot.png")
  # save the plot as png file if 'bSaveFigures' is set to be TRUE
  bSaveFigures <- TRUE
  if(bSaveFigures==T){
    ggsave(plt08,file=fgNm08,
           width=210,height=297*1.4,
           #width=210,height=297,
           #width=297,height=210,
           units="mm",dpi=300)
  }
  
  
  #______________________________________________________________________________
  # Linear regression 01 - start
  #______________________________________________________________________________
  # --------- regression ----------------------
  # make a list of species abbreviations 
  ls.spcAbbr <-  unique(df_e06.sb$speciesabbr)
  no.ofspcs <- length(ls.spcAbbr) 
  #assign to a different variable
  iter <- no.ofspcs
  #make a variable that defines number of columns for a matrix
  vars = 4
  #prepare an empty matrix w enough rows
  mtx_plts2 <- matrix(ncol=vars, nrow=no.ofspcs)
  # set a number for the growing number of elements in the list
  k <- 1
  #iterate over species
  for (spclat in ls.spcAbbr){
    print(spclat)
  #}
  # subset data frame based on species abbreviation
  sbs_csvs03 <- df_e06.sb[df_e06.sb$speciesabbr==spclat,]
  #Exclude too low and to high values, as they are outliers for the ddPCR machine
  #calculate the covariance
  cov_sbs02 <- cov(log10(sbs_csvs03$mcp_qPCR+1), 
                   log10(sbs_csvs03$mcp_ddPCR+1))
  #calculate the correlation
  cor_sbs02 <- cor(-log10(sbs_csvs03$mcp_qPCR+1), 
                   log10(sbs_csvs03$mcp_ddPCR+1))*100
  rcor_sbs02 <- round(cor_sbs02, 3)
  rcor_sbs02[is.na(rcor_sbs02)] <- 0
  #estimate a linear model 
  lin.md01 <-  lm(log10(mcp_qPCR+1)~log10(mcp_ddPCR+1),sbs_csvs03)
  #get the slope to calculate the efficiency
  slo1 <- lin.md01$coefficients[2]
  slo2 <- as.numeric(as.character(slo1))
  # get intercept
  intc1 <- lin.md01$coefficients[1]
  intc2 <- as.numeric(as.character(intc1))
  intc3 <- round(intc2,3)
  # get slope
  slo3 <- round(slo2,3)
  RSqdRn <- rcor_sbs02
  # add to matrix
  mtx_plts2[k,] <- c((as.character(spclat)),
                     intc3,slo3,RSqdRn)
  # add to increasing number
  k <- k+1
  # end iteration over species
  }
#make the matrix a data frame
df_lm.mod02 <- as.data.frame(mtx_plts2)
# change column names
colnames(df_lm.mod02) <- c("spcAbbr",
                           "intercpt",
                           "slope",
                           "RSqdRn")
write.csv(df_lm.mod02, 
          file=paste0(wd00_01,"/table_w_linear_model_regression_results_v06.csv"))

#______________________________________________________________________________
# Linear regression 01 - end
#______________________________________________________________________________
#______________________________________________________________________________
# Linear regression 02 -  start 2nd way of getting regression data
#_______________________________________________________________________________
library(purrr)
library(dplyr)
library(tidyr)
# get regression data for both MST and STD samples
df <- df_e06.sb %>% 
  dplyr::rename(ddPCR=mcp_ddPCR, qPCR=mcp_qPCR) %>%
  dplyr::mutate(ddPCR=ifelse(ddPCR==0,NA,ddPCR)) %>%
  dplyr::mutate(qPCR=ifelse(qPCR==0,NA,qPCR))

df_regr <- df %>%
  dplyr::filter(!is.na(ddPCR),!is.na(qPCR)) %>%
  split(.$speciesabbr) %>%
  purrr::map(~ lm(log10(qPCR) ~ log10(ddPCR), data = .))

df_fit <- df_regr %>%
  purrr::map_dfr(broom::tidy,.id="speciesabbr") 

df_fit <- df_fit %>%
  dplyr::select(speciesabbr,term,estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate)

df_stats <- df_regr %>%
  purrr::map(summary) %>%
  purrr::map_dfr(broom::glance,.id="speciesabbr") %>%
  dplyr::select(speciesabbr,`r.squared`,`adj.r.squared`,`p.value`)

df_fit <- df_fit %>%
  dplyr::left_join(df_stats,by="speciesabbr")
#Write out the table
# https://stackoverflow.com/questions/7303322/apply-function-to-each-column-in-a-data-frame-observing-each-columns-existing-da
# apply function to all columns if they are numeric, otherwise return the column content
mtx_fit2 <- sapply( df_fit, function(x) if("numeric" %in% class(x) ) { 
  round(as.numeric(as.character(x)),3)
} else { (x) } )
# make the matrix a data frame
df_fit2 <- as.data.frame(mtx_fit2)
df_fit2$p.value <-  df_fit$p.value

write.csv(df_fit2, file=paste0(wd00_01,
                               "/table_w_linear_model_regression_results_v07.csv"))
# get regression data for only STD samples
df <- df_e06.sb %>% 
  dplyr::rename(ddPCR=mcp_ddPCR, qPCR=mcp_qPCR) %>%
  dplyr::mutate(ddPCR=ifelse(ddPCR==0,NA,ddPCR)) %>%
  dplyr::mutate(qPCR=ifelse(qPCR==0,NA,qPCR))

df_regr <- df %>%
  dplyr::filter(!is.na(ddPCR),!is.na(qPCR)) %>%
  split(.$speciesabbr) %>%
  purrr::map(~ lm(log10(qPCR) ~ log10(ddPCR), data = .))

df_fit <- df_regr %>%
  purrr::map_dfr(broom::tidy,.id="speciesabbr") 

df_fit <- df_fit %>%
  dplyr::select(speciesabbr,term,estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate)

df_stats <- df_regr %>%
  purrr::map(summary) %>%
  purrr::map_dfr(broom::glance,.id="speciesabbr") %>%
  dplyr::select(speciesabbr,`r.squared`,`adj.r.squared`,`p.value`)

df_fit <- df_fit %>%
  dplyr::left_join(df_stats,by="speciesabbr")

# https://stackoverflow.com/questions/7303322/apply-function-to-each-column-in-a-data-frame-observing-each-columns-existing-da
# apply function to all columns if they are numeric, otherwise return the column content
mtx_fit2 <- sapply( df_fit, function(x) if("numeric" %in% class(x) ) { 
  round(as.numeric(as.character(x)),3)
} else { (x) } )
# make the matrix a data frame
df_fit2 <- as.data.frame(mtx_fit2)
df_fit2$p.value <-  df_fit$p.value

#Write out the table
write.csv(df_fit2, file=paste0(wd00_01,"/table_w_linear_model_regression_results_v08.csv"))

#_______________________________________________________________________________
# linear regression 02 - end-  2nd way of getting regression data
#_______________________________________________________________________________

# #_______________________________________________________________________________

# #_______________________________________________________________________________
# get the phyla and class for the species
df_e06.sb$ClNm <-  df_dtc_asss$Class[match(df_e06.sb$speciesabbr,df_dtc_asss$AbbrvNm)]
df_e06.sb$PhNm <-  df_dtc_asss$Phyla[match(df_e06.sb$speciesabbr,df_dtc_asss$AbbrvNm)]


# read in file with MST locations for MST samples  
inf05 <- "table03_combined_MST_samples_locations.csv"
wd00.1_inf05 <- paste0(wd00.1,"/",inf05 )
library(xlsx)
df_MSTsmpls <- read.csv2(wd00.1_inf05,sep=";")
# limit to only comprise samples from 2021
df_MSTsmpls <- df_MSTsmpls[grepl("2021-",df_MSTsmpls$Dato_inds),]
# check if all are FALSE - in this case there are no duplicates
all(duplicated(df_MSTsmpls$U_Pr_Nr))
# use the dplyr library to only select a couple of the columns
library(dplyr)
df_MSTsmpls02 <- df_MSTsmpls %>% dplyr::select(U_Pr_Nr,
                              Dato_inds,
                              lokalitet_vanda,
                              lok_pos_lat,
                              lok_pos_lon,
                              Lok_omr01)

# copy a column , and substitute inside the values of this column, 
# so that there is a common column with the same identifiers
# to use for 'left_join'
df_e06.sb$U_Pr_Nr <- gsub("(^[0-9]{+})-","\\1",df_e06.sb$smplNm)
# use left_join to get the    'lok_pos_lat', 'lok_pos_lon'
# added on to the data frame that holds the ddPCR and the qPCR data
# store this in a new data frame
df_e07 <- df_e06.sb %>% dplyr::left_join(df_MSTsmpls02,
                               by="U_Pr_Nr")

# copy the columns with lat and lon to have them appearing in the
# data frame with shorter names, to have shorter names, when they are to
# be called in the plotting
df_e07$lat <- as.numeric(df_e07$lok_pos_lat)
df_e07$lon <- as.numeric(df_e07$lok_pos_lon)
# Try making the bubble plot as found on the website
#

#https://r-graph-gallery.com/330-bubble-map-with-ggplot2.html
# Libraries
library(ggplot2)
library(dplyr)
# Get the world polygon and extract the country
#get giscoR package
if(!require(giscoR)){
  #https://ropengov.github.io/giscoR/
  #remotes::install_github("rOpenGov/giscoR")
  install.packages("giscoR",repos = c("https://ropengov.r-universe.dev", 
                                      "https://cloud.r-project.org")
  )
}  
library(giscoR)
DK_map <- gisco_get_countries(country = c("Denmark",
                                          "Sweden",
                                          "Norway",
                                          "Germany",
                                          "Poland"
                                          ), resolution = 1)
# Get a data frame with longitude, latitude, and size of bubbles (a bubble = a city)
library(maps)
# PLot the collection points
pl07 <- ggplot() +
  geom_sf(data = DK_map, fill = "grey34", alpha = 0.3) +
  geom_point(data = df_e07, aes(x = lon, 
                                y = lat)) +
  #theme_void() +
  xlim(4, 17) +
  ylim(54, 59)
# only keep the columns required for the mapping
df_e08 <- df_e07 %>% dplyr::select(speciesabbr,
                                   smplNm,
                                   U_Pr_Nr,
                                   mcp_ddPCR,
                                   mcp_qPCR,
                                   GnNm,
                                   SpNm,
                                   ClNm,
                                   PhNm,
                                   speciesNm,
                                   lokalitet_vanda,
                                   lon,
                                   lat,
                                   Lok_omr01,
                                   Dato_inds)
# re arrange the data frame to  a long format with the copy count of'
# of DNA molecules
df_e08 <- df_e08 %>% tidyr::pivot_longer(.,
                    c(mcp_ddPCR,mcp_qPCR), 
                    names_to = c('typCnt','mch'), 
                    values_to = "molcnt",
                    #names_pattern = '(.*)\\.(t\\d+)') %>%
                    names_pattern = '(.*)_(.*)') %>%
  dplyr::group_by(speciesabbr , smplNm)
# Get teh collection dates and get the day the month and the year
df_e08$Dato_inds.dd <- as.numeric(gsub("([0-9]{4})-([0-9]{2})-([0-9]{2})",
                                       "\\3",df_e08$Dato_inds))
df_e08$Dato_inds.mm <- as.numeric(gsub("([0-9]{4})-([0-9]{2})-([0-9]{2})",
                                       "\\2",df_e08$Dato_inds))
df_e08$Dato_inds.yy <- as.numeric(gsub("([0-9]{4})-([0-9]{2})-([0-9]{2})",
                                       "\\1",df_e08$Dato_inds))
# evaluate on the month and assign a season category,
# to an NA filled column
df_e08$ssn.smpl <- NA

# df_e08$ssn.smpl[(df_e08$Dato_inds.mm>7)] <- "1st_season_Jan_to_Jun"
# df_e08$ssn.smpl[(df_e08$Dato_inds.mm<7)] <- "2nd_season_Jul_to_Nov"

df_e08$ssn.smpl[(df_e08$Dato_inds.mm>7)] <- "1st_season"
df_e08$ssn.smpl[(df_e08$Dato_inds.mm<7)] <- "2nd_season"
library(RColorBrewer)
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)      # Show all color palettes

# evaluate on the season, if there is no season, then the
# row is not needed for the plot that is to be prepared
df_e09 <- df_e08[!is.na(df_e08$ssn.smpl),]
# make a column that holds the combined information on
# the season sampled and the machine
df_e09$mch.ssn <- paste0(df_e09$mch," : ",df_e09$ssn.smpl)
# make a column that has log10 evaluations of the molecular count
df_e09$eDN.lv1.l10 <- log10(df_e09$molcnt+1)
# calcaluate how to skew sampling locations in a grid
uspcAb  <- unique(df_e09$speciesabbr)
nuspcAb <- length(uspcAb)
npa <- round(sqrt(nuspcAb),0)
ppa <- ceiling(nuspcAb/npa)
sppa <- seq(1,ppa,1)
snpa <- seq(1,npa,1)

snpa <- seq(-npa/2,(npa/2)-1,1)
sppa <- seq(-ppa/2,(ppa/2)-1,1)
npb <- npa*ppa
snpb <- seq(1,npb,1)
lon.m <- rep(sppa,npa)
lat.m <- rep(snpa,ppa)
df_coord.m <- as.data.frame(cbind(snpb,lon.m,lat.m))
df_coord.m <- df_coord.m[1:nuspcAb,]
df_coord.m <- cbind(df_coord.m,uspcAb)
df_coord.m$lon.m <- df_coord.m$lon.m*0.16
df_coord.m$lat.m <- df_coord.m$lat.m*0.08
df_coord.m$speciesabbr <- df_coord.m$uspcAb
df_e10 <- df_e09 %>% dplyr::left_join(df_coord.m,
                                      by="speciesabbr")
# adjust the lon and lat
df_e10$lon.m <- df_e10$lon.m+df_e10$lon
df_e10$lat.m <- df_e10$lat.m+df_e10$lat
# no need to include the not monitored species
df_e10 <- df_e10[(df_e10$molcnt!=0),]
library(ggrepel)
library(ggplot2)

# make plot with sampled locations
pl07 <- ggplot() +
  geom_sf(data = DK_map, fill = "grey71", alpha = 1.0) +
  geom_point(data = df_e10, aes(x = lon, 
                                y = lat),
             
             size=11,
             fill=alpha("white",alpha = 0),
             colour='black',
             #alpha = 0.0,
             stroke=0.2,
             shape=22) +
  
  geom_point(data = dplyr::arrange(df_e10,across(eDN.lv1.l10,desc)),
             aes(x = lon.m, 
                 y = lat.m,
                 colour=mch,
                 fill = speciesabbr),
             size=1.2,
             alpha = 0.7,
             stroke=0.2,
             shape=21) +
  theme_bw() +
  #https://stackoverflow.com/questions/62009919/facet-title-alignment-using-facet-wrap-in-ggplot2
  theme(strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
        panel.background = element_rect(fill = "white",
                                        color = "white")) +
  
  
  scale_fill_viridis_d(option = 'viridis') +
  #scale_color_manual(values=c("springgreen4", "#E69F00")) +
  scale_color_manual(values=c("black", "black")) +
  labs(color='machine') +
  labs(size='log10(copy/uL)') +
  labs(fill='species') +
  facet_grid(mch ~ ssn.smpl) +
  # see : https://www.datanovia.com/en/blog/how-to-change-ggplot-facet-labels/
  theme(strip.background = element_rect(
    color="white", 
    fill="white", size=1.5, linetype="solid")
  ) +
  theme(
    strip.text.x = element_text(
      size = 10, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, color = "black", face = "bold"
    )
  ) +
  #facet_wrap(mch ~ ssn.smpl, ncol = 2) +
  guides(colour="none") +
  guides(size="none") +
  #theme_void() +
  xlim(7.6, 16.8) +
  ylim(54.2, 58) +
  xlab("longitude") + ylab("latitude")

pl07
fgNm08 <- paste0(wd00_01,"/Fig13_v01_map_ddPCR_qPCR.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- TRUE
if(bSaveFigures==T){
  ggsave(pl07,file=fgNm08,
         #width=210,height=297*1.4,
         width=210,height=297*0.4,
         #width=297,height=210,
         units="mm",dpi=300)
}


# make plot with sampled locations
pl07 <- ggplot() +
  geom_sf(data = DK_map, fill = "grey71", alpha = 1.0) +
  geom_point(data = df_e10, aes(x = lon, 
                                y = lat),
             
             size=11,
             fill=alpha("white",alpha = 0),
             colour='black',
             #alpha = 0.0,
             stroke=0.2,
             shape=22) +
  
  geom_point(data = dplyr::arrange(df_e10,across(eDN.lv1.l10,desc)),
             aes(x = lon.m, 
                 y = lat.m,
                 colour=mch,
                 size=eDN.lv1.l10,
                 fill = speciesabbr),
             #size=1.2,
             alpha = 0.7,
             stroke=0.2,
             shape=21) +
  theme_bw() +
  #https://stackoverflow.com/questions/62009919/facet-title-alignment-using-facet-wrap-in-ggplot2
  theme(strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
        panel.background = element_rect(fill = "white",
                                        color = "white")) +
  
  
  scale_fill_viridis_d(option = 'viridis') +
  #scale_color_manual(values=c("springgreen4", "#E69F00")) +
  scale_color_manual(values=c("black", "black")) +
  labs(color='machine') +
  labs(size='log10(copy/uL)') +
  labs(fill='species') +
  facet_grid(mch ~ ssn.smpl) +
  # see : https://www.datanovia.com/en/blog/how-to-change-ggplot-facet-labels/
  theme(strip.background = element_rect(
    color="white", 
    fill="white", size=1.5, linetype="solid")
  ) +
  theme(
    strip.text.x = element_text(
      size = 10, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, color = "black", face = "bold"
    )
  ) +
  #facet_wrap(mch ~ ssn.smpl, ncol = 2) +
  guides(colour="none") +
  #guides(size="none") +
  #theme_void() +
  xlim(7.6, 16.8) +
  ylim(54.2, 58) +
  xlab("longitude") + ylab("latitude")

pl07
fgNm08 <- paste0(wd00_01,"/Fig13_v02_map_ddPCR_qPCR.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- TRUE
if(bSaveFigures==T){
  ggsave(pl07,file=fgNm08,
         #width=210,height=297*1.4,
         width=210,height=297*0.5,
         #width=297,height=210,
         units="mm",dpi=300)
}
#_______________________________________________________________________
# get the unique individual phyla names
# to be able to iterate over them
uPhNm <- unique(df_e10$PhNm)
nPhNm <- length(uPhNm)
sPhNm <- seq(1,nPhNm,1)
# iterate over numbers for phyla names
for (phNM in sPhNm)
{
  # subset the data frame
  sbsPhnm <- uPhNm[phNM]
  df_e11 <- df_e10[(df_e10$PhNm==sbsPhnm),]
  # make plot with sampled locations
  
  pl07 <- ggplot() +
    geom_sf(data = DK_map, fill = "grey71", alpha = 1.0) +
    geom_point(data = df_e11, aes(x = lon,
                                  y = lat),

               size=11,
               fill=alpha("white",alpha = 0),
               colour='black',
               #alpha = 0.0,
               stroke=0.2,
               shape=22) +
    
    geom_point(data = dplyr::arrange(df_e11,across(eDN.lv1.l10,desc)),
               aes(x = lon.m, 
                   y = lat.m,
                   colour=mch,
                   size=eDN.lv1.l10,
                   fill = speciesabbr),
               #size=1.2,
               alpha = 0.7,
               stroke=0.2,
               shape=21) +
    theme_bw() +
    #https://stackoverflow.com/questions/62009919/facet-title-alignment-using-facet-wrap-in-ggplot2
    theme(strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
          panel.background = element_rect(fill = "white",
                                          color = "white")) +
    
    
    scale_fill_viridis_d(option = 'viridis') +
    #scale_color_manual(values=c("springgreen4", "#E69F00")) +
    scale_color_manual(values=c("black", "black")) +
    labs(color='machine') +
    labs(size='log10(copy/uL)') +
    labs(fill='species') +
    facet_grid(mch ~ ssn.smpl) +
    # see : https://www.datanovia.com/en/blog/how-to-change-ggplot-facet-labels/
    theme(strip.background = element_rect(
      color="white", 
      fill="white", size=1.5, linetype="solid")
    ) +
    theme(
      strip.text.x = element_text(
        size = 10, color = "black", face = "bold"
      ),
      strip.text.y = element_text(
        size = 10, color = "black", face = "bold"
      )
    ) +
    #facet_wrap(mch ~ ssn.smpl, ncol = 2) +
    guides(colour="none") +
    #guides(size="none") +
    #theme_void() +
    xlim(7.6, 16.8) +
    ylim(54.2, 58) +
    xlab("longitude") + ylab("latitude")
  
  pl07
  fgNm08 <- paste0(wd00_01,"/Fig14_v0",phNM,"_map_ddPCR_qPCR_",sbsPhnm,".png")
  # save the plot as png file if 'bSaveFigures' is set to be TRUE
  bSaveFigures <- TRUE
  if(bSaveFigures==T){
    ggsave(pl07,file=fgNm08,
           #width=210,height=297*1.4,
           width=210,height=297*0.5,
           #width=297,height=210,
           units="mm",dpi=300)
  }
}
#
# https://www.urbandemographics.org/post/figures-map-layers-r/