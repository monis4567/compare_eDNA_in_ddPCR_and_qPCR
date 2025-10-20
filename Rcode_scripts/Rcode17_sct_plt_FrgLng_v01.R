#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
#remove everything in the working environment, without a warning!!
#rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)


library("dplyr") # Load dplyr package
library("plyr") # Load plyr package
library("readr") # Load readr package


#library("rnaturalearthhires")
#define working directory
#define input file  directory
wd01 <- "ddpcr_resultater"
wd01 <- "data_ddpcr_runs"
wd00 <- getwd()
wd_d <- "data"
#setwd(wd00)
getwd()
# paste together working directory and input file  directory
wd00_01 <- paste0(wd00,"/",wd_d,"/",wd01)
#grep for the '.csv' files in the directory, and place in a list
lst_fcsv <- list.files(wd00_01, full.names = T)[grep("\\.csv",list.files(wd00_01))]
# get all files that has 
lst_fcsv <- lst_fcsv[grep("results_from_",lst_fcsv)]
lst_fcsv <- lst_fcsv[grep("all_assays",lst_fcsv)]
# change the working directory to where the csv files are
#setwd(wd00_01)
# follow this example
# https://statisticsglobe.com/merge-csv-files-in-r
# To read in all files
df_csvs01 <- lst_fcsv %>%
  lapply(read_csv) %>% # Store all files in list
  bind_rows # Combine data sets into one data set 
# change the working directory again
#setwd(wd00)
# define directory for original merged qPCR textreports
wdin01<- "MONIS6_2021_data/output02_merged_txtfiles_from_mxpro_for_MONIS6"
# define input flie to read in
inf01 <- "outfile02_merged_mxpro_csvfls_MONIS6.csv"

wdin01_inf01 <- paste(wd00,wd_d,wdin01,inf01,sep="/")
# read in delimited file with all merged qPCR results
df_qPCRm01 <- read_delim(wdin01_inf01,delim=";")

# read in an xlsx file with all primer assays listed 
# for each species and the abbreviations
wd00.1 =paste0(wd00,"/data/MONIS6_2021_data")
inf04 <- "list_of_specific_assays_MONIS6.xlsx"
wd00.1_inf04 <- paste0(wd00.1,"/",inf04 )
df_dtc_asss <- openxlsx::read.xlsx(wd00.1_inf04,1)
# paste together working directory and output file  directory
wdout <- paste0(wd00,"/output11_compare_all_ddpcr_and_qpcr")
wd05 <- "output04_stdcrv_plots_and_tables_from_Rcode_for_MST2017_2023_samples"
# define path for previously prepared table
#for standard curve efficiencies
pth_and_fl2 <- paste(wd00,"/",wd05,"/table05_2_MONIS6_eDNA_std_crv_efficiencies.csv", sep="")
# read in the ';' separated file , read the header row as column names
df_stdef <- read.csv(pth_and_fl2,sep=";",header=T)
# Use the 'spec.lat' column to get the genus name and species name
# by splitting in to two columns using the ' ' space as delimiter
df_stdef <- df_stdef %>%
  tidyr::separate(spec.lat, into = c("genus", "species"), sep = " ")
# get the first 3 letters from the string in the 'genus' column
df_stdef$gen3 <- substr(df_stdef$genus,1,3)
# get the first 3 letters from the string in the 'species' column
df_stdef$spc3 <- substr(df_stdef$species,1,3)
# paste together the first 3 letters of genus and species to get a species abbreviation
df_stdef$speciesabbr <- paste0(df_stdef$gen3,df_stdef$spc3)


#___________________________________
# calculate LOD and LOQ for qPCR data
#___________________________________
# 
df_qPCRm01$snm_pltn <- paste0(df_qPCRm01$speciesabbr,"_",df_qPCRm01$plateno)
#get the LOD for each qPCR run
#use the function aggregate to get the minimum value for a group
lodtable1 <- aggregate(df_qPCRm01[, "Quantitycopies"], list(df_qPCRm01$snm_pltn, df_qPCRm01$WellType), min)
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

#add an empty column with just NAs
df_qPCRm01[,"eDNA_eval"] <- NA
#replace in the empty column, the order is important, as you otherwise will end up with the last evaluations
df_qPCRm01$eDNA_eval[df_qPCRm01$QuanCp02>=df_qPCRm01$LOQ] <- "aboveLOQ"
df_qPCRm01$eDNA_eval[df_qPCRm01$QuanCp02<df_qPCRm01$LOQ] <- "AbLOD_BeLOQ"
df_qPCRm01$eDNA_eval[df_qPCRm01$QuanCp02<df_qPCRm01$LOD] <- "belowLOD"
df_qPCRm01$eDNA_eval[df_qPCRm01$QuanCp02==0] <- "NoCt"
# evaluate the eDNA levels
# check if the eDNA level is zero
df_qPCRm01[,"eDNA_NC"] <- 0
df_qPCRm01$eDNA_NC[df_qPCRm01$QuanCp02==0] <- 1
# check if the eDNA level is  below LOD and also not zero
df_qPCRm01[,"eDNA_bLD"] <- 0
df_qPCRm01$eDNA_bLD[df_qPCRm01$QuanCp02<df_qPCRm01$LOD & df_qPCRm01$QuanCp02!=0] <- 1

# check if the eDNA level is above LOD and also below LOQ
df_qPCRm01[,"eDNA_aLDbLQ"] <- 0
df_qPCRm01$eDNA_aLDbLQ[df_qPCRm01$QuanCp02<df_qPCRm01$LOQ & df_qPCRm01$QuanCp02>=df_qPCRm01$LOD] <- 1
# check if the eDNA level is above LOQ
df_qPCRm01[,"eDNA_aLQ"] <- 0
df_qPCRm01$eDNA_aLQ[df_qPCRm01$QuanCp02>=df_qPCRm01$LOQ] <- 1

library(dplyr)
library(plyr)

# Count evaluations per sample per species
# https://stackoverflow.com/questions/28090119/summing-all-columns-by-group
df_qPCRm02 <- df_qPCRm01 %>%  plyr::summarise(smpltp,speciesabbr,plateno,
                                              eDNA_NC,eDNA_bLD,eDNA_aLDbLQ,eDNA_aLQ) %>% 
  dplyr::group_by(smpltp,speciesabbr, plateno) %>%
  dplyr::summarise_each(list(sum))
# add a color column for evaluations of eDNA levels inferred by qPCR 
# per sample per species
df_qPCRm02$cl.f.evl <- "NA"
# evaluate for each sample and species
df_qPCRm02$cl.f.evl[df_qPCRm02$eDNA_NC>=3 ] <- "white"
df_qPCRm02$cl.f.evl[df_qPCRm02$eDNA_bLD>=1 ] <- "yellow"
df_qPCRm02$cl.f.evl[df_qPCRm02$eDNA_aLDbLQ>=1 ] <- "orange"
df_qPCRm02$cl.f.evl[df_qPCRm02$eDNA_aLQ>=1 ] <- "red"
df_qPCRm02$cl.f.evl[df_qPCRm02$eDNA_aLQ>=3 ] <- "black"
#paste together species abbreviation and sample name
df_qPCRm02$spA.smplN <- paste0(df_qPCRm02$speciesabbr,".",df_qPCRm02$smpltp)
#___________________________________
df_csvs01$spcAbbr <- df_csvs01$speciesabbr
# substitute wrong species abbreviation
df_csvs01$spcAbbr <- gsub("Myaara","Myaare",df_csvs01$spcAbbr)
df_csvs01$smplTp <- substring(df_csvs01$smplNm,1,3)
# only retain MST samples
df_csvs02 <- df_csvs01[df_csvs01$smplTp=="MST",]
# copy column
df_csvs02$f2 <- df_csvs02$sp.f2
df_csvs02$mcp_ddPCR_lod <- df_csvs02$ddlod
# evaluate whether the copy levels are above the LOD for the platform used
df_csvs02$evl1.qP <- df_csvs02$mcp_qPCR>(df_csvs02$qlod*df_csvs02$f2)
df_csvs02$evl1.dP <- df_csvs02$mcp_ddPCR>df_csvs02$mcp_ddPCR_lod
#paste together species abbreviation and sample name
df_csvs02$spA.smplN <- paste0(df_csvs02$spcAbbr,".",df_csvs02$smplNm)

# define column names to keep
keep.colm1 <- c("smplNm",
                "spcAbbr",
                "mcp_qPCR",
                "mcp_ddPCR")
# only retain specified columns
df_csvs03 <- df_csvs02[keep.colm1]
# define column names to keep
keep.colm2 <- c("smplNm",
                "spcAbbr",
                "evl1.qP","evl1.dP")
# only retain specified columns
df_csvs03.2 <- df_csvs02[keep.colm2]
# exclude any spcAbbr that might be NAs
df_csvs03<-df_csvs03[!is.na(df_csvs03$spcAbbr),]
# use tidyr to re arrange the data frame
df_csvs04 <- df_csvs03 %>%
  tidyr::pivot_wider(names_from="spcAbbr",
                     values_from=c("mcp_qPCR","mcp_ddPCR"))
# transform bolean opreators to ones and zeros
df_csvs03.2$evl1.qP <- df_csvs03.2$evl1.qP*1
df_csvs03.2$evl1.dP <- df_csvs03.2$evl1.dP*1
df_csvs03.3 <- df_csvs03.2 %>% dplyr::count(spcAbbr,evl1.qP, evl1.dP)
library(ggplot2)
# count the columns
ncl4 <- ncol(df_csvs04)
#View(df_csvs04)
#rearrange the data frame
df_csvs05 <- df_csvs04 %>% 
  tidyr::pivot_longer(c(2:ncl4), names_to="sample") 

# exclude the failed and empty ddPCR attempts
df_csvs05 <- df_csvs05[!grepl("Neomel",df_csvs05$sample),]
df_csvs05 <- df_csvs05[!grepl("Oncgor",df_csvs05$sample),]
df_csvs05 <- df_csvs05[!grepl("Rhihar",df_csvs05$sample),]

# modify sample name
df_csvs05$sample <- gsub("mcp_","",df_csvs05$sample)
# swap sample ID around
df_csvs05$sample <- gsub("^(.*)_(.*)$","\\2_\\1",df_csvs05$sample)
df_csvs05$machine <- gsub("^(.*)_(.*)$","\\2",df_csvs05$sample)
# make a list with colors for the sampling categories
clfH <- rep(c("grey54","white"),length(unique(df_csvs05$smplNm)))

# get species and make breaks
arter <-sort(unique(df_csvs05$sample))
arter_breaks <- arter[seq(1,length(arter),2)]
arter_labels <- stringr::str_remove(arter_breaks,"_ddPCR")
# make range of intervals to use for horizontal lines in the plot
# split string by delimiter to get species name
df_csvs05$spcNm <- sapply(strsplit(df_csvs05$sample,"_"), "[[", 1)
# dfv01 <- data.frame(do.call('rbind',
#             strsplit(as.character(df_csvs05$value),
#                      ',',fixed=TRUE)))
# vl1 <-  as.numeric(gsub("c\\(","",dfv01$X1))
# vl2 <-  gsub("\\)","",dfv01$X2)
# vl2 <-  as.numeric(gsub(" ","",vl2))
# df_csvs05$vl1 <- vl1
# 
# df_csvs05$vl2 <- vl2
# value.vec <- df_csvs05$value
# define columns to keep
ctke    <- c("smplNm",
             "sample",
             "machine",
             "spcNm")
#
df_csvs05 <- df_csvs05[ctke]

nospc <- length(unique(df_csvs05$spcNm))
nosmp <- length(unique(df_csvs05$smplNm))
# make range of intervals to use for horizontal lines in the plot
y_values_hline <- seq(2.5, (2*(nospc)+0.5), 2)
# make range of intervals to use for vertical lines in the plot
x_values_hline <- seq(3.5, (3*(nosmp)+0.5), 3)


# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# df_csvs05 <- lapply(df_csvs05, as.numeric)
# # replace zeroes with NAs to avoid getting a point in the bubble chart
# df_csvs05$vl1[(df_csvs05$vl1==0)] <- NA
# df_csvs05$vl2[(df_csvs05$vl2==0)] <- NA
# 
# #use 'tidyr::pivot_longer' to get a long version of the data frame
# nvlcolm <- colnames(df_csvs05)[grepl("vl",colnames(df_csvs05))]
# df_csvs05.1 <- df_csvs05 %>% tidyr::pivot_longer(cols=c(nvlcolm),
#                                                  names_to = "vlc", 
#                                                  values_to="vln")
# # paste together to get a column with both species name and sample name
# df_csvs05.1$spA.smplN <- paste0(df_csvs05.1$spcNm,".",df_csvs05.1$smplNm)
# # match to get colour for evaluation of eDNA level at sample 
# # location per species
# df_csvs05.1$cl.f.evl <- df_qPCRm02$cl.f.evl[match(df_csvs05.1$spA.smplN,df_qPCRm02$spA.smplN)]
# # add white for NAs
# df_csvs05.1$cl.f.evl[is.na(df_csvs05.1$cl.f.evl)] <- "white"
# # copy column
# df_csvs05.1$eval.c1 <- df_csvs05.1$cl.f.evl
# # copy data frame
# df_csvs05 <- df_csvs05.1
# #
# gplot4 <- ggplot(data=df_csvs05.1, aes(smplNm, sample, 
#                                        color = vlc, size = log10(vln))) + 
#   theme_classic() +
#   geom_point(aes(colour=machine)) +
#   scale_size(range = c(0.1, 10), name="log10(molecules DNA/uL)") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   # use section below to get rectangel backgrounds - however
#   # this will require a second fill scale - and I was not able to 
#   # get this second color scale working
#   # geom_tile(aes(height=Inf,fill=as.factor(smplNm),y=0.005),
#   #           position=position_nudge(x=-0.05), show.legend=F) +
#   # scale_fill_manual(values=alpha(clfH,0.04)) +
#   # # # 
#   scale_color_manual(values=c("springgreen4", "#E69F00")) +
#   geom_hline(yintercept = y_values_hline, colour = "gray64") +
#   labs(color='machine')+
#   scale_y_discrete(limits=rev) +
#   xlab("MST sample") + ylab("species")
# # see the ggplot object
# #gplot4
# 
# # set parameter to check whether to save plot or not
# bSaveFigures <- T
# # make a filename with path
# fnm02 <- paste0(wd00_01,"/Fig05_v02_compare_all_ddPCR_and_qPCR.png")
# #p
# if(bSaveFigures==T){
#   ggsave(gplot4,file=fnm02,
#          #width=297*0.6,height=210*1.2,
#          height=297*0.6,width=210,
#          units="mm",dpi=300)
# }
# #
# df_csvs05$eDN.lv1 <- df_csvs05$vln
# df_csvs05$eDN.lv1[is.na(df_csvs05$eDN.lv1)] <- 0
# df_csvs05$eDN.lv1 <- as.numeric(df_csvs05$eDN.lv1)
# # copy columns
# df_csvs01$spcAbbr <- df_csvs01$speciesabbr
# df_csvs01$f2 <- df_csvs01$sp.f2
# df_csvs01$mcp_ddPCR_lod <- df_csvs01$ddlod
# # get unique values for the species abbreviations
# # for the LOD calculated per ddPCR analysis
# df_qlod.ddP01 <- df_csvs01 %>% dplyr::group_by(spcAbbr) %>%
#   distinct(spcAbbr, mcp_ddPCR_lod,f2,qlod)
# # use match to transfer the qlod and factor for difference between qPCR
# # and ddPCR
# df_csvs05$f2 <- df_qlod.ddP01$f2[match(df_csvs05$spcNm,
#                                        df_qlod.ddP01$spcAbbr)]
# df_csvs05$dd.qlod <- df_qlod.ddP01$qlod[match(df_csvs05$spcNm,
#                                               df_qlod.ddP01$spcAbbr)]
# df_csvs05$mcp_ddPCR_lod <- df_qlod.ddP01$mcp_ddPCR_lod[match(
#   df_csvs05$spcNm,df_qlod.ddP01$spcAbbr)]
# # calculate loq10 levels of eDNA levels. #
# # Add 1 to ensure log10 to zero is not infinite
# df_csvs05$eDN.lv1.l10 <- log10(df_csvs05$eDN.lv1+1)
# # now remove all lgo10 to 1, as these are zero, and are not needed on the plot
# df_csvs05$eDN.lv1.l10[df_csvs05$eDN.lv1.l10==0] <- NA
# # subset the dataframe to only comprise ddPCR findings
# df_dP06 <- df_csvs05[(df_csvs05$machine=="ddPCR"),]
# # copy the column with evaluation for colors
# df_dP06$eval.c2 <- df_dP06$eval.c1
# 
# # evaluate a color coding for eDNA levels on ddPCR in relation to the qlod 
# # and mcp_ddPCR_lod
# df_dP06$eval.c2[df_dP06$eDN.lv1==0] <- "white"
# df_dP06$eval.c2[df_dP06$eDN.lv1!=0 & df_dP06$eDN.lv1<df_dP06$dd.qlod] <- "yellow"
# df_dP06$eval.c2[df_dP06$eDN.lv1>df_dP06$dd.qlod] <- "black"
# # append back to combined data frame
# df_csvs05$eval.c1[(df_csvs05$machine=="ddPCR")] <- df_dP06$eval.c2
# # make a path for wrting a temporary file
# wd00_01_tmpf05 <- paste0(wd00_01,"/tmp_df_csvs05.csv")
# # write a temporary file for making a reprex
# write.csv(df_csvs05, wd00_01_tmpf05)
# # copy the evaluation column
# df_csvs05$eval.c4 <- df_csvs05$eval.c1
# #modify the evaluation category names
# df_csvs05$eval.c1[df_csvs05$eval.c1=="white"] <- "no Cq" # no eDNA
# df_csvs05$eval.c1[df_csvs05$eval.c1=="yellow"] <- "bLD" # below LOD
# df_csvs05$eval.c1[df_csvs05$eval.c1=="orange"] <- "aLDbLQ" # above LOD below LOQ
# df_csvs05$eval.c1[df_csvs05$eval.c1=="red"] <- "1aLQ" # 1 above LOQ
# df_csvs05$eval.c1[df_csvs05$eval.c1=="black"] <- "aaLQ" # all above LOQ
# # count the number of assays tested
# no_of_as_Tst <- length(unique(df_csvs05$spcNm))
# # ensure the plot library is loaded
# library(ggplot2)
# # begin the plot
# plt_01a <- ggplot(df_csvs05, aes(smplNm, 
#                                  sample, 
#                                  size = eDN.lv1.l10,
#                                  color = machine, 
#                                  fill = eval.c1)) +
#   theme_classic() +
#   geom_point(shape = 21, stroke = 1.2) +
#   geom_hline(yintercept = y_values_hline, colour = alpha("gray64",0.6)) +
#   geom_vline(xintercept = x_values_hline, colour = alpha("gray64",0.6)) +
#   scale_color_manual(values=alpha(c("springgreen4", "#E69F00"),0.9)) +
#   scale_fill_manual(values=alpha(
#     c("red","black","orange","yellow","white")
#     ,0.5),name="evaluation") +
#   scale_size(range = c(0.001, 10), 
#              name="log10(molecules\n DNA/uL)") +
#   scale_y_discrete(name="assay", 
#                    breaks=arter_breaks, labels=arter_labels) +
#   scale_x_discrete(name="watersample") +
#   coord_flip() +
#   labs(color='machine') +
#   theme(axis.text.x = element_text(angle = 90, 
#                                    vjust = 0.5, hjust=1)) 
# # adjust the size of the header categories on the legend key text
# plt_01a  <- plt_01a + theme(legend.title=element_text(size=10))
# #make a file to save to
# fnm02 <- paste0(wd00_01,"/Fig05_v03_compare_all_ddPCR_and_qPCR.png")
# #evaluate if the plot should be saved
# if(bSaveFigures==T){
#   ggsave(plt_01a,file=fnm02,
#          #width=297*0.6,height=210*1.2,
#          height=297*0.6,width=210,
#          units="mm",dpi=300)
# }
# 
# setwd(wd00_01)
# # --------- regression ----------------------
# # make a list of species abbreviations 
# ls.spcAbbr <-  unique(df_csvs01$spcAbbr)
# no.ofspcs <- length(ls.spcAbbr) 
# #assign to a different variable
# iter <- no.ofspcs
# #make a variable that defines number of columns for a matrix
# vars = 4
# #prepare an empty matrix w enough rows
# mtx_plts2 <- matrix(ncol=vars, nrow=no.ofspcs)
# # set a number for the growing number of elements in the list
# k <- 1
# #iterate over species
# for (spclat in ls.spcAbbr){
#   print(spclat)
#   # subset data frame based on species abbreviation
#   sbs_csvs03 <- df_csvs01[df_csvs01$spcAbbr==spclat,]
#   #Exclude too low and to high values, as they are outliers for the ddPCR machine
#   #calculate the covariance
#   cov_sbs02 <- cov(log10(sbs_csvs03$mcp_qPCR+1), 
#                    log10(sbs_csvs03$mcp_ddPCR+1))
#   #calculate the correlation
#   cor_sbs02 <- cor(-log10(sbs_csvs03$mcp_qPCR+1), 
#                    log10(sbs_csvs03$mcp_ddPCR+1))*100
#   rcor_sbs02 <- round(cor_sbs02, 3)
#   rcor_sbs02[is.na(rcor_sbs02)] <- 0
#   #estimate a linear model 
#   lin.md01 <-  lm(log10(mcp_qPCR+1)~log10(mcp_ddPCR+1),sbs_csvs03)
#   #get the slope to calculate the efficiency
#   slo1 <- lin.md01$coefficients[2]
#   slo2 <- as.numeric(as.character(slo1))
#   # get intercept
#   intc1 <- lin.md01$coefficients[1]
#   intc2 <- as.numeric(as.character(intc1))
#   intc3 <- round(intc2,3)
#   # get slope
#   slo3 <- round(slo2,3)
#   RSqdRn <- rcor_sbs02
#   # add to matrix
#   mtx_plts2[k,] <- c((as.character(spclat)),
#                      intc3,slo3,RSqdRn)
#   # add to increasing number
#   k <- k+1
#   # end iteration over species
# }
# #make the matrix a data frame
# df_lm.mod02 <- as.data.frame(mtx_plts2)
# # change column names
# colnames(df_lm.mod02) <- c("spcAbbr",
#                            "intercpt",
#                            "slope",
#                            "RSqdRn")
# write.csv(df_lm.mod02, file=paste0(wd00_01,"/table_w_linear_model_regression_results_v02.csv"))
# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#_______________________________________________________________________________
# start 2nd way of getting regression data
#_______________________________________________________________________________
library(purrr)
library(dplyr)
library(tidyr)
# get regression data for both MST and STD samples
df <- df_csvs01 %>% 
  dplyr::rename(ddPCR=mcp_ddPCR, qPCR=mcp_qPCR) %>%
  dplyr::mutate(ddPCR=ifelse(ddPCR==0,NA,ddPCR)) %>%
  dplyr::mutate(qPCR=ifelse(qPCR==0,NA,qPCR))

df_regr <- df %>%
  dplyr::filter(!is.na(ddPCR),!is.na(qPCR)) %>%
  split(.$spcAbbr) %>%
  purrr::map(~ lm(log10(qPCR) ~ log10(ddPCR), data = .))

df_fit <- df_regr %>%
  purrr::map_dfr(broom::tidy,.id="spcAbbr") 

df_fit <- df_fit %>%
  dplyr::select(spcAbbr,term,estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate)

df_stats <- df_regr %>%
  purrr::map(summary) %>%
  purrr::map_dfr(broom::glance,.id="spcAbbr") %>%
  dplyr::select(spcAbbr,`r.squared`,`adj.r.squared`,`p.value`)

df_fit <- df_fit %>%
  dplyr::left_join(df_stats,by="spcAbbr")
#Write out the table
# https://stackoverflow.com/questions/7303322/apply-function-to-each-column-in-a-data-frame-observing-each-columns-existing-da
# apply function to all columns if they are numeric, otherwise return the column content
mtx_fit2 <- sapply( df_fit, function(x) if("numeric" %in% class(x) ) { 
  round(as.numeric(as.character(x)),3)
} else { (x) } )
# make the matrix a data frame
df_fit2 <- as.data.frame(mtx_fit2)
df_fit2$p.value <-  df_fit$p.value

#write.csv(df_fit2, file=paste0(wd00_01,"/table_w_linear_model_regression_results_v03.csv"))
# get regression data for only STD samples
df <- df_csvs03 %>% 
  dplyr::rename(ddPCR=mcp_ddPCR, qPCR=mcp_qPCR) %>%
  dplyr::mutate(ddPCR=ifelse(ddPCR==0,NA,ddPCR)) %>%
  dplyr::mutate(qPCR=ifelse(qPCR==0,NA,qPCR))

df_regr <- df %>%
  dplyr::filter(!is.na(ddPCR),!is.na(qPCR)) %>%
  split(.$spcAbbr) %>%
  purrr::map(~ lm(log10(qPCR) ~ log10(ddPCR), data = .))

df_fit <- df_regr %>%
  purrr::map_dfr(broom::tidy,.id="spcAbbr") 

df_fit <- df_fit %>%
  dplyr::select(spcAbbr,term,estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate)

df_stats <- df_regr %>%
  purrr::map(summary) %>%
  purrr::map_dfr(broom::glance,.id="spcAbbr") %>%
  dplyr::select(spcAbbr,`r.squared`,`adj.r.squared`,`p.value`)

df_fit <- df_fit %>%
  dplyr::left_join(df_stats,by="spcAbbr")

# https://stackoverflow.com/questions/7303322/apply-function-to-each-column-in-a-data-frame-observing-each-columns-existing-da
# apply function to all columns if they are numeric, otherwise return the column content
mtx_fit2 <- sapply( df_fit, function(x) if("numeric" %in% class(x) ) { 
  round(as.numeric(as.character(x)),3)
} else { (x) } )
# make the matrix a data frame
df_fit2 <- as.data.frame(mtx_fit2)
df_fit2$p.value <-  df_fit$p.value

#Write out the table
#write.csv(df_fit2, file=paste0(wd00_01,"/table_w_linear_model_regression_results_v04.csv"))

#_______________________________________________________________________________
# end 2nd way of getting regression data
#_______________________________________________________________________________
# exclude  negative controls
df_csvs01.1 <- df_csvs01[(df_csvs01$smplTp!="NTC"),]
df_csvs01.1 <- df_csvs01.1[(df_csvs01.1$smplTp!="NEK"),]
# use dplyr to rename columns
df <- df_csvs01.1 %>%
  dplyr::rename(ddPCR=mcp_ddPCR,
                qPCR=mcp_qPCR)

# begin ggplot
p <- ggplot(df) +
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=smplTp,
                 shape=smplTp,
                 colour=smplTp), size=2) +
  # geom_errorbar(aes(ymin = sdmn_qPCR,ymax = sdmx_qPCR)) + 
  # geom_errorbarh(aes(ymin = sdmn_ddPCR,ymax = sdmx_ddPCR)) +
  geom_smooth(aes(x=ddPCR, y=qPCR),
              method="lm",
              alpha=0.4,
              fill="#CCCCCC",
              #fill="green",
              #colour="#FF0000",
              colour="blue",
              #colour=smplTp,
              size=0.5) +
  scale_alpha_ordinal() +
  #facet_wrap(.~spcAbbr, ncol=2, scales="free") +
  facet_wrap(.~spcAbbr, ncol=4) +
  theme_minimal() +
  scale_fill_manual(values=c("red","pink")) +
  scale_color_manual(values=c(rep("black",2))) +
  scale_shape_manual(values=c(rep(21,2))) +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth  =0.5)) +
  xlab("ddPCR, molecules DNA/uL") + ylab("qPCR, molecules DNA/uL") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  labs(color='sample') +
  labs(shape='sample') +
  labs(fill='sample') +
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
#p

# make a 2nd version of the linear model plot
df <- df_csvs01.1 %>% 
  dplyr::rename(ddPCR=mcp_ddPCR, qPCR=mcp_qPCR) %>%
  mutate(ddPCR=ifelse(ddPCR==0,NA,ddPCR)) %>%
  mutate(qPCR=ifelse(qPCR==0,NA,qPCR))

list_species <- sort(unique(df$spcAbbr))
df$spcAbbr <- factor(df$spcAbbr,levels=list_species)
# same axes limits across facets
p2 <- ggplot(df) +
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=smplTp,
                 shape=smplTp,
                 colour=smplTp), size=2,show.legend=T) +
  geom_smooth(aes(x=ddPCR, y=qPCR),
              method="lm",
              alpha=0.4,
              fill="#CCCCCC",
              #colour="#FF0000",
              colour="blue",
              size=0.5) +
  facet_wrap(.~spcAbbr, ncol=4, drop=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth =0.5)) +
  xlab("ddPCR, molecules DNA/uL") + ylab("qPCR, molecules DNA/uL") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  labs(color='sample') +
  labs(shape='sample') +
  labs(fill='sample') +
  
  scale_fill_manual(values=c("red", "pink")) +
  scale_color_manual(values=c(rep("black",2))) +
  scale_shape_manual(values=c(21,24)) +
  
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
# make the legend position be on top of the plot
p2 <- p2 + theme(legend.position = "top")

# get the first 3 characters of the string
df_dtc_asss$gn.abbr <- substr(df_dtc_asss$Genus,1,3)
df_dtc_asss$sp.abbr <- substr(df_dtc_asss$Species,1,3)
# paste together the two strings with no space
df_dtc_asss$speciesabbr <- paste0(df_dtc_asss$gn.abbr,df_dtc_asss$sp.abbr)


# exclude if there is no 'Fprimseq'
# if it is NA
df_das <- df_dtc_asss[!is.na(df_dtc_asss$Fprimseq),]
# substitute in the colunm names so that 'seq' and 'Seq'
# are the same, by being replaced with 'seq'
colnames(df_das) <- gsub("Seq","seq",colnames(df_das))

# use left_join to add 'df_das' to the  'df_e06.sb'
# data frame
df2 <- left_join(df, 
                 df_das, by=c("speciesabbr"="speciesabbr"))

# exclude rows of the target fragment length is above 150
#df2 <- df2[df2$TargetFraglngt <= 150,]
# exclude rows when 'smplTp' is MST
#df2 <- df2[!df2$smplTp=="MST",]

# same axes limits across facets
p2 <- ggplot(df2) +
  # plot the error bars first so they are below the points
  geom_errorbar(aes(x=ddPCR, y=qPCR,ymin =sdmn_qPCR  ,ymax =sdmx_qPCR  ),
                orientation = "x") + 
  geom_errorbar(aes(x=ddPCR, y=qPCR,xmin =sdmn_ddPCR,xmax = sdmx_ddPCR),
                orientation = "y") +
  # next plot the points
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=TargetFraglngt,
                 shape=smplTp,
                 colour=smplTp), size=7.1,show.legend=T) +
  geom_smooth(aes(x=ddPCR, y=qPCR),
              method="lm",
              alpha=0.4,
              fill="#CCCCCC",
              #colour="#FF0000",
              colour="blue",
              size=0.5) +
  #facet_wrap(.~spcAbbr, ncol=4, drop=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth =0.5),
        axis.text = element_text(size = 18)) +  # Adjust tick label size here
  xlab("ddPCR, (molecules DNA/uL)") + ylab("qPCR, (molecules DNA/uL)") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
labs(color='sample') +
  labs(shape='sample') +
  labs(fill='fragment\nlength (bp)') +
  #scale_fill_manual(values=c("red", "pink")) +
  scale_color_manual(values=c(rep("black",2))) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_viridis_c(option = "F", direction = -1, alpha = 0.7) +  # Set alpha for fill here
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
# make the legend position be on top of the plot
#p2 <- p2 + theme(legend.position = "top")
# make the legend position be on bottom of the plot
p2 <- p2 + theme(legend.position = "right")

fgN <- paste0(wdout,"/Fig14_v02_comb_scatterplot.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- T
if(bSaveFigures==T){
  ggsave(p2,file=fgN,
         width=210*1.2,height=297*0.8,
         #width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}

# use left_join to add 'df_das' to the  'df_e06.sb'
# data frame
df2 <- left_join(df, 
                 df_das, by=c("speciesabbr"="speciesabbr"))

# exclude rows if both the column 'qPCR' and 'ddPCR' are NA
df2 <- df2[!is.na(df2$qPCR) | !is.na(df2$ddPCR),]
# exclude rows if either the column 'qPCR' and 'ddPCR' are NA
df2 <- df2[!is.na(df2$qPCR) & !is.na(df2$ddPCR),]

# exclude rows of the target fragment length is above 150
#df2 <- df2[df2$TargetFraglngt <= 150,]
# exclude rows when 'smplTp' is MST
#df2 <- df2[!df2$smplTp=="MST",]
# subset the data frame to only include rows where 'speciesabbr' 
# is "Mnelei" or "Myaare" or "Psever or "Psefar"
df2 <- df2[df2$speciesabbr %in% c("Mnelei",
                                  "Myaare",
                                  #"Promin",
                                  #"Karmik",
                                  "Psever",
                                  "Psefar"),]
#View(df2)
# Make a plot
p2 <- ggplot(df2) +
  # plot the error bars first so they are below the points
  geom_errorbar(aes(x=ddPCR, y=qPCR,ymin =sdmn_qPCR  ,ymax =sdmx_qPCR  ),
                orientation = "x") + 
  geom_errorbar(aes(x=ddPCR, y=qPCR,xmin =sdmn_ddPCR,xmax = sdmx_ddPCR),
                orientation = "y") +
  # next plot the points
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=TargetFraglngt,
                 shape=smplTp,
                 colour=smplTp), size=7.1,show.legend=T) +
  # geom_smooth(aes(x=ddPCR, y=qPCR),
  #             method="lm",
  #             alpha=0.4,
  #             fill="#CCCCCC",
  #             #colour="#FF0000",
  #             colour="blue",
  #             size=0.5) +
  #facet_wrap(.~spcAbbr, ncol=4, drop=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth =0.5),
  axis.text = element_text(size = 18)) +  # Adjust tick label size here
  xlab("ddPCR, (molecules DNA/uL)") + ylab("qPCR, (molecules DNA/uL)") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
labs(color='sample') +
  labs(shape='sample') +
  labs(fill='fragment\nlength (bp)') +
  
  #scale_fill_manual(values=c("red", "pink")) +
  scale_color_manual(values=c(rep("black",2))) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_viridis_c(option = "E", direction = -1, alpha = 0.7) +  # Set alpha for fill here
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
# make the legend position be on bottom of the plot
p2 <- p2 + theme(legend.position = "right")
p2
fgN <- paste0(wdout,"/Fig14_v03_comb_scatterplot.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- T
if(bSaveFigures==T){
  ggsave(p2,file=fgN,
         width=210*1.2,height=297*0.8,
         #width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}

# exclude rows when 'smplTp' is not MST
#df2 <- df2[df2$smplTp=="MST",]

# Make a plot
p2 <- ggplot(df2) +
  # plot the error bars first so they are below the points
  geom_errorbar(aes(x=ddPCR, y=qPCR,ymin =sdmn_qPCR  ,ymax =sdmx_qPCR  ),
                orientation = "x") + 
  geom_errorbar(aes(x=ddPCR, y=qPCR,xmin =sdmn_ddPCR,xmax = sdmx_ddPCR),
                orientation = "y") +
  # next plot the points
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=TargetFraglngt, shape=smplTp),
             size=7.1,show.legend=T) +
  # geom_smooth(aes(x=ddPCR, y=qPCR),
  #             method="lm",
  #             alpha=0.4,
  #             fill="#CCCCCC",
  #             #colour="#FF0000",
  #             colour="blue",
  #             size=0.5) +
  #facet_wrap(.~spcAbbr, ncol=4, drop=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth =0.5),
        axis.text = element_text(size = 18)) +  # Adjust tick label size here
  xlab("ddPCR, (molecules DNA/uL)") + ylab("qPCR, (molecules DNA/uL)") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
labs(color='sample') +
  labs(shape='sample') +
  labs(fill='fragment\nlength (bp)') +
  
  #scale_fill_manual(values=c("red", "pink")) +
  #scale_color_manual(values=c(rep("black",2))) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_viridis_c(option = "B", direction = -1, alpha = 0.7) +  # Set alpha for fill here
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
# make the legend position be on bottom of the plot
p2 <- p2 + theme(legend.position = "right")
Fig14_v04 <- p2
fgN <- paste0(wdout,"/Fig14_v04_comb_sctplt_wtout_std.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- T
if(bSaveFigures==T){
  ggsave(p2,file=fgN,
         width=210*1.2,height=297*0.8,
         #width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}

#########################################################################
#
# Calculate Tm, ΔG, and GC content for primers and probes
#
#########################################################################

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998)
nn_params <- data.frame(
  dinuc = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"),
  dH     = c(-7.9, -8.4, -7.8, -7.2, -8.5, -8.0, -10.6, -7.8, -8.2, -9.8, -8.0, -8.4, -7.2, -8.2, -8.5, -7.9),
  dS     = c(-22.2, -22.4, -21.0, -20.4, -22.7, -19.9, -27.2, -21.0, -22.2, -24.4, -19.9, -22.4, -21.3, -22.2, -22.7, -22.2)
)

# 1. Function to Calculate Tm (Melting Temperature)
calculate_tm <- function(dna_seq, Na_conc = 40, primer_conc = 50e-5) {
  seq <- toupper(dna_seq)
  n <- nchar(seq)
  if (n == 0) return(NA)
  
  # Initialize
  dH <- 0
  dS <- 0
  
  # Calculate ΔH and ΔS
  for (i in 1:(n-1)) {
    dinuc <- paste(substr(seq, i, i), substr(seq, i+1, i+1), sep = "")
    params <- nn_params[nn_params$dinuc == dinuc, ]
    if (nrow(params) > 0) {
      dH <- dH + params$dH[1]  # Take the first match
      dS <- dS + params$dS[1]
    }
  }
  
  # Add initiation and symmetry corrections
  dH <- dH + 0.2  # initiation
  dS <- dS - 5.7  # initiation
  
  # Salt correction
  R <- 1.987  # Gas constant in cal/mol·K
  salt_corr <- 16.6 * log10(Na_conc / 1000)
  
  # Calculate Tm
  tm <- (1000 * dH) / (dS + R * log(primer_conc)) - 273.15 + salt_corr
  
  return(tm)
}

# 2. Function to Calculate ΔG (Gibbs Free Energy) at a Given Temperature
calculate_deltaG <- function(dna_seq, Na_conc = 40, primer_conc = 50e-5, temp_C = 60) {
  seq <- toupper(dna_seq)
  n <- nchar(seq)
  if (n == 0) return(NA)
  
  # Initialize
  dH <- 0
  dS <- 0
  
  # Calculate ΔH and ΔS
  for (i in 1:(n-1)) {
    dinuc <- paste(substr(seq, i, i), substr(seq, i+1, i+1), sep = "")
    params <- nn_params[nn_params$dinuc == dinuc, ]
    if (nrow(params) > 0) {
      dH <- dH + params$dH[1]  # Take the first match
      dS <- dS + params$dS[1]
    }
  }
  
  # Add initiation and symmetry corrections
  dH <- dH + 0.2  # initiation
  dS <- dS - 5.7  # initiation
  
  # Convert temperature to Kelvin
  temp_K <- temp_C + 273.15
  
  # Calculate ΔG in kcal/mol
  deltaG <- dH - (temp_K / 1000) * dS
  
  return(deltaG)
}

# 3. Function to Calculate GC Content
calculate_gc_content <- function(dna_seq) {
  seq <- toupper(dna_seq)
  n <- nchar(seq)
  if (n == 0) return(NA)
  
  gc <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
gc_content <- (gc / n) * 100

return(gc_content)
}

# Assign row names
row.names(df_das) <- 1:nrow(df_das)

# Use the functions to calculate Tm, ΔG, and GC content for primers and probes
# Define a list to hold calculated variables
results_list <- list()

# Iterate over each row of the data frame
for (i in 1:nrow(df_das)) {
  # Access specific row using the current iteration number
  current_row <- df_das[i, ]
  
  # Calculate variables
  FprimTm <- calculate_tm(df_das$Fprimseq[i])
  RprimTm <- calculate_tm(df_das$Rprimseq[i])
  ProbeTm <- calculate_tm(df_das$Probeseq[i])
  
  # Calculate ΔG at 60°C
  FprimDeltaG <- calculate_deltaG(df_das$Fprimseq[i], temp_C = 60)
  RprimDeltaG <- calculate_deltaG(df_das$Rprimseq[i], temp_C = 60)
  ProbeDeltaG <- calculate_deltaG(df_das$Probeseq[i], temp_C = 60)
  
  # Calculate GC content
  FprimGC <- calculate_gc_content(df_das$Fprimseq[i])
  RprimGC <- calculate_gc_content(df_das$Rprimseq[i])
  ProbeGC <- calculate_gc_content(df_das$Probeseq[i])
  
  # Store results in the list
  results_list[[i]] <- data.frame(
    FprimTm = FprimTm,
    RprimTm = RprimTm,
    ProbeTm = ProbeTm,
    FprimDeltaG = FprimDeltaG,
    RprimDeltaG = RprimDeltaG,
    ProbeDeltaG = ProbeDeltaG,
    FprimGC = FprimGC,
    RprimGC = RprimGC,
    ProbeGC = ProbeGC
  )
}

# Combine the results list into a data frame
results_df <- do.call(rbind, results_list)

# Merge the results back to original data frame
df_das <- cbind(df_das, results_df)

# View the updated data frame
# View(df_das)
# Inspired by this paper:
# https://link.springer.com/article/10.1007/s10529-022-03295-2
# make a scatter plot w ggplot 
# GC content of the probe on the x-axis
# and the DeltaG of the probe on the y-axis
p3 <- ggplot(df_das) +
  # make poitns for the probe
  geom_point(aes(x=ProbeGC, y=ProbeDeltaG),
             size=3, color="blue", fill="lightblue", shape=21) +
  # make poitns for the R primer
  geom_point(aes(x=RprimGC, y=RprimDeltaG),
             size=3, color="red", fill="pink", shape=22) +
  # make poitns for the F primer
  geom_point(aes(x=FprimGC, y=FprimDeltaG),
             size=3, color="green", fill="lightgreen", shape=24) +
  theme_minimal() +
  xlab("GC Content of Probe (%)") +
  ylab(expression(Delta* G~"(kcal/mol)")) +
  ggtitle("Scatter Plot of GC Content vs. ΔG of oligo") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -20, linetype="dashed", color = "red") +
  geom_vline(xintercept = 50, linetype="dashed", color = "red") +
  annotate("text", x = max(df_das$ProbeGC, na.rm=TRUE)*0.8, 
           y = -18, label = "ΔG = -20 kcal/mol", color = "red") +
  annotate("text", x = 52, 
           y = max(df_das$ProbeDeltaG, na.rm=TRUE)*0.8, 
           label = "GC Content = 50%", color = "red")
p3



# use left_join to add 'df_das' to the  'df_e06.sb'
# data frame
df3 <- left_join(df, 
                 df_das, by=c("speciesabbr"="speciesabbr"))

# exclude rows if both the column 'qPCR' and 'ddPCR' are NA
df3 <- df3[!is.na(df3$qPCR) | !is.na(df3$ddPCR),]
# exclude rows if either the column 'qPCR' and 'ddPCR' are NA
df3 <- df3[!is.na(df3$qPCR) & !is.na(df3$ddPCR),]

# subset the data frame to only include rows where 'speciesabbr' 
# is "Mnelei" or "Myaare" or "Psever or "Psefar"
df3 <- df3[df3$speciesabbr %in% c("Mnelei",
                                  "Myaare",
                                  #"Promin",
                                  #"Karmik",
                                  "Psever",
                                  "Psefar"),]
#View(df2)
# Make a plot
p2 <- ggplot(df3) +
  # plot the error bars first so they are below the points
  geom_errorbar(aes(x=ddPCR, y=qPCR,ymin =sdmn_qPCR  ,ymax =sdmx_qPCR  ),
                orientation = "x") + 
  geom_errorbar(aes(x=ddPCR, y=qPCR,xmin =sdmn_ddPCR,xmax = sdmx_ddPCR),
                orientation = "y") +
  
  # next plot the points
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=ProbeDeltaG, shape = smplTp),
              size=7.1,show.legend=T) +
  # geom_smooth(aes(x=ddPCR, y=qPCR),
  #             method="lm",
  #             alpha=0.4,
  #             fill="#CCCCCC",
  #             #colour="#FF0000",
  #             colour="blue",
  #             size=0.5) +
  #facet_wrap(.~spcAbbr, ncol=4, drop=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth =0.5),
        axis.text = element_text(size = 18)) +  # Adjust tick label size here
  xlab("ddPCR, (molecules DNA/uL)") + ylab("qPCR, (molecules DNA/uL)") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
labs(color='sample') +
  labs(shape='sample') +
  #labs(fill=expression("Probe"~Delta* G)) +
  labs(fill="Probe\nΔG") +
  
  
  #scale_fill_manual(values=c("red", "pink")) +
  #scale_color_manual(values=c(rep("black",2))) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_viridis_c(option = "D", direction = -1, alpha = 0.7,
                       labels = scales::number_format(accuracy = 0.1)) +  # Set alpha for fill here
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
# make the legend position be on bottom of the plot
p2 <- p2 + theme(legend.position = "right")
Fig14_v05 <- p2
fgN <- paste0(wdout,"/Fig14_v05_comb_sctplt_wtout_std.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- T
if(bSaveFigures==T){
  ggsave(p2,file=fgN,
         width=210*1.2,height=297*0.8,
         #width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}




# use left_join to add 'df_das' to the  'df'
# data frame
df3 <- left_join(df, 
                 df_das, by=c("speciesabbr"="speciesabbr"))


# use left_join to add 'df_stdef' to the  'df3'
# data frame
df4 <- left_join(df3, 
                 df_stdef, by=c("speciesabbr"="speciesabbr"))

#View(df4)
# exclude rows if both the column 'qPCR' and 'ddPCR' are NA
df4 <- df4[!is.na(df4$qPCR) | !is.na(df4$ddPCR),]
# exclude rows if either the column 'qPCR' and 'ddPCR' are NA
df4 <- df4[!is.na(df4$qPCR) & !is.na(df4$ddPCR),]



# # subset the data frame to only include rows where 'speciesabbr'
# # is "Mnelei" or "Myaare" or "Psever or "Psefar"
# df4 <- df4[df4$speciesabbr %in% c("Mnelei",
#                                   "Myaare",
#                                   #"Promin",
#                                   #"Karmik",
#                                   "Psever",
#                                   "Psefar"),]
# add all deltaG values for each oligo
df4$oligoDeltaG <- df4$FprimDeltaG + df4$RprimDeltaG + df4$ProbeDeltaG
# Make a plot
p2 <- ggplot(df4) +
  # plot the error bars first so they are below the points
  geom_errorbar(aes(x=ddPCR, y=qPCR,ymin =sdmn_qPCR  ,ymax =sdmx_qPCR  ),
                orientation = "x") + 
  geom_errorbar(aes(x=ddPCR, y=qPCR,xmin =sdmn_ddPCR,xmax = sdmx_ddPCR),
                orientation = "y") +
  
  # next plot the points
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=oligoDeltaG),
             shape=21, size=7.1,show.legend=T) +
  # geom_smooth(aes(x=ddPCR, y=qPCR),
  #             method="lm",
  #             alpha=0.4,
  #             fill="#CCCCCC",
  #             #colour="#FF0000",
  #             colour="blue",
  #             size=0.5) +
  #facet_wrap(.~spcAbbr, ncol=4, drop=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth =0.5),
        axis.text = element_text(size = 18)) +  # Adjust tick label size here
  xlab("ddPCR, (molecules DNA/uL)") + ylab("qPCR, (molecules DNA/uL)") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
labs(color='sample') +
  labs(shape='sample') +
  #labs(fill=expression("oligo "~Delta* " G")) +
  labs(fill="oligo\nΔG") +
  
  
  #scale_fill_manual(values=c("red", "pink")) +
  #scale_color_manual(values=c(rep("black",2))) +
  #scale_shape_manual(values=c(21,24)) +
  scale_fill_viridis_c(option = "A", direction = -1, alpha = 0.7) +  # Set alpha for fill here
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
# make the legend position be on bottom of the plot
p2 <- p2 + theme(legend.position = "right")
p2
fgN <- paste0(wdout,"/Fig14_v06_comb_sctplt_w_deltaG.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- T
if(bSaveFigures==T){
  ggsave(p2,file=fgN,
         width=210*1.2,height=297*0.8,
         #width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}



# use left_join to join the dataframes  'df_das' and 'df_stdef'
# by the column 'speciesabbr'
df_das2 <- left_join(df_das,
                    df_stdef, by=c("speciesabbr"="speciesabbr"))
# make sure the column 'rEffic' is numeric
df_das2$rEffic <- as.numeric(df_das2$rEffic)
# subset the 'df_das2' with 'dplyr' to only comprise the columns
# speciesabbr, FprimDeltaG, RprimDeltaG, ProbeDeltaG, rEffic, 
# TargetFraglngt, FprimGC, RprimGC, ProbeGC,
# FprimTm, RprimTm, ProbeTm, Fprimseq, RprimSeq, ProbeSeq
# FprimLngt, RprimLngt, ProbeLngt, Fprimseq, RprimSeq, ProbeSeq
df_das3 <- df_das2 %>%
  dplyr::select(speciesabbr, FprimDeltaG, RprimDeltaG, ProbeDeltaG,
         rEffic, TargetFraglngt,
         FprimGC, RprimGC, ProbeGC,
         FprimTm, RprimTm, ProbeTm,
         #FprimLngt, RprimLngt, ProbeLngt,
         Fprimseq, Rprimseq, Probeseq)


# make the data frame longer with pivot longer to have a column
# that starts with 'oligo'  for the columns 
# that currently start with 'Fprim', 'Rprim', and 'Probe'
df_das_long <- df_das3 %>%
  pivot_longer(cols = c(FprimDeltaG, RprimDeltaG, ProbeDeltaG,
                        FprimGC, RprimGC, ProbeGC,
                        FprimTm, RprimTm, ProbeTm,
                        Fprimseq, Rprimseq, Probeseq),
               names_to = c("oligo", ".value"),
               names_pattern = "(Fprim|Rprim|Probe)(.*)")
#View the long data frame
df_das_long$oligo <- as.factor(df_das_long$oligo)

# Inspired by this paper:
# https://link.springer.com/article/10.1007/s10529-022-03295-2
# make a scatter plot w ggplot 
# GC content of the probe on the x-axis
# and the DeltaG of the probe on the y-axis
p3 <- ggplot() +
  # make poitns for the probe
  geom_point(data=df_das_long,
             aes(x=GC, y=DeltaG,
                 color=oligo, 
                 fill=rEffic),
             stroke=0.5,
             size=3, 
             shape=21) +
  theme_minimal() +
  xlab("GC Content of Probe (%)") +
  ylab(expression(Delta* G~"(kcal/mol)")) +
  ggtitle("Scatter Plot of GC Content vs. ΔG of oligo") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -20, linetype="dashed", color = "red") +
  geom_vline(xintercept = 50, linetype="dashed", color = "red") +
  scale_fill_viridis_c(option = "D", direction = -1, alpha = 0.7) +  # Set alpha for fill here
  labs(color='oligo') +
  labs(fill='rEffic') +
  scale_color_manual(values=c("red", "blue", "orange")) +
  annotate("text", x = max(df_das$ProbeGC, na.rm=TRUE)*0.8, 
           y = -18, label = "ΔG = -20 kcal/mol", color = "red") +
  annotate("text", x = 52, 
           y = max(df_das$ProbeDeltaG, na.rm=TRUE)*0.8, 
           label = "GC Content = 50%", color = "red")
p3

fgN <- paste0(wdout,"/Fig14_v08_GC_Gibbs_energ.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- T
if(bSaveFigures==T){
  ggsave(p3,file=fgN,
         width=210*1.2,height=297*0.8,
         #width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}






# use left_join to add 'df_das' to the  'df'
# data frame
df3 <- left_join(df, 
                 df_das2, by=c("speciesabbr"="speciesabbr"))

unique(df$smplTp)

#View(df)
# exclude rows if both the column 'qPCR' and 'ddPCR' are NA
df4 <- df4[!is.na(df4$qPCR) | !is.na(df4$ddPCR),]
# exclude rows if either the column 'qPCR' and 'ddPCR' are NA
df4 <- df4[!is.na(df4$qPCR) & !is.na(df4$ddPCR),]



# subset the data frame to only include rows where 'speciesabbr'
# is "Mnelei" or "Myaare" or "Psever or "Psefar"
df4 <- df4[df4$speciesabbr %in% c("Mnelei",
                                  "Myaare",
                                  #"Promin",
                                  #"Karmik",
                                  "Psever",
                                  "Psefar"),]


min(df4$TargetFraglngt)
max(df4$TargetFraglngt)

min(df4$ProbeDeltaG)
max(df4$ProbeDeltaG)
# Make a plot
p2 <- ggplot(df4) +
  # plot the error bars first so they are below the points
  geom_errorbar(aes(x=ddPCR, y=qPCR,ymin =sdmn_qPCR  ,ymax =sdmx_qPCR  ),
                orientation = "x") + 
  geom_errorbar(aes(x=ddPCR, y=qPCR,xmin =sdmn_ddPCR,xmax = sdmx_ddPCR),
                orientation = "y") +
  # next plot the points
  geom_point(aes(x=ddPCR, y=qPCR,
                 fill=rEffic, shape=smplTp),
             size=7.1,show.legend=T) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA,
                                    colour="#999999",
                                    linewidth =0.5),
        axis.text = element_text(size = 18)) +  # Adjust tick label size here
  xlab("ddPCR, (molecules DNA/uL)") + ylab("qPCR, (molecules DNA/uL)") +
  theme(
    axis.text.x = element_text(size = 17),    # Adjust font size of x-axis tick labels
    axis.text.y = element_text(size = 17),    # Adjust font size of y-axis tick labels
    axis.title.x = element_text(size = 16),  # Adjust size as needed
    axis.title.y = element_text(size = 16),   # Adjust size as needed
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  labs(color='sample') +
  labs(shape='sample') +
  labs(fill='rEffic') +
  #scale_fill_viridis_c(option = "D", direction = -1, alpha = 0.7) +  # Set alpha for fill here
  scale_fill_gradient2(
    low = adjustcolor("darkblue", alpha.f = 0.7),   # Semi-transparent blue for values below 100
    mid = adjustcolor("white", alpha.f = 0.7),  # Semi-transparent white for value 100
    high = adjustcolor("red", alpha.f = 0.7),    # Semi-transparent red for values above 100
    midpoint = 100,
    na.value = "grey",
    guide = "colorbar"   # Optional: adds a color bar
  ) +
  scale_shape_manual(values = c(21,24)) +  # Manually setting shapes
  annotation_logticks(sides="bl") +
  geom_abline(intercept = 0, slope = 1, 
              color="darkorchid3",linetype = "dashed") +
  scale_y_log10(labels = label_log(),limits=c(1e-2, 1e5)) +
  scale_x_log10(labels = label_log(),limits=c(1e-2, 1e5))
# make the legend position be on bottom of the plot
p2 <- p2 + theme(legend.position = "right")
Fig14_v07 <- p2
fgN <- paste0(wdout,"/Fig14_v07_comb_sctplt_w_rEffic.png")
# save the plot as png file if 'bSaveFigures' is set to be TRUE
bSaveFigures <- T
if(bSaveFigures==T){
  ggsave(p2,file=fgN,
         width=210*1.2,height=297*0.8,
         #width=210,height=297,
         #width=297,height=210,
         units="mm",dpi=300)
}
library(cowplot)
# Combine plots into one column
combined_plot <- 
  cowplot::plot_grid( Fig14_v04, 
                      Fig14_v05, 
                      Fig14_v07,
                      nrow=3,ncol = 1,
                     labels = c('(a)', '(b)','(c)'),
                     label_size = 26)
 
fgN <- paste0(wdout,"/Fig14_v09_combine_3_plts.png")
# Save the combined plot
ggsave(
  fgN, 
  plot = combined_plot, 
  width = 14.8, # half the width of A4 in cm
  height = 29.7, # full height of A4 in cm
  dpi = 300
)
