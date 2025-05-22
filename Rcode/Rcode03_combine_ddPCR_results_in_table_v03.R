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


library("dplyr") # Load dplyr package
library("plyr") # Load plyr package
library("readr") # Load readr package
# install package
if(!require(rnaturalearth)){
  install.packages("rnaturalearth")
  install.packages("devtools")
  #https://github.com/ropensci/rnaturalearthhires
  remotes::install_github("ropensci/rnaturalearthhires")
  install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")
}
#library("rnaturalearthhires")
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
# get all files that has 
lst_fcsv <- lst_fcsv[grep("results_from_",lst_fcsv)]
lst_fcsv <- lst_fcsv[grep("all_assays",lst_fcsv)]
lst_fcsv <- lst_fcsv[grep("02",lst_fcsv)]
# change the working directory to where the csv files are
setwd(wd00_01)
# follow this example
# https://statisticsglobe.com/merge-csv-files-in-r
# To read in all files
df_csvs01 <- lst_fcsv %>%
  lapply(read_csv) %>% # Store all files in list
  bind_rows # Combine data sets into one data set 
# Check if any rows are duplicated for the columns "speciesabbr","smplNm"
which(duplicated(df_csvs01[c("speciesabbr","smplNm")]))
# change the working directory again
setwd(wd00)
# define directory for original merged qPCR textreports
wdin01<- paste0(wd00,"/MONIS6_2021_data/output02_merged_txtfiles_from_mxpro_for_MONIS6")

# define input flie to read in
inf01 <- "outfile02_merged_mxpro_csvfls_MONIS6.csv"
wdin01_inf01 <- paste0(wdin01,"/",inf01)
# read in delimited file with all merged qPCR results
df_qPCRm01 <- read_delim(wdin01_inf01,delim=";")


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
#             strsplit(as.character(df_csvs05$value),',',fixed=TRUE)))
# colnames(dfv01)
# vl1 <-  as.numeric(gsub("c\\(","",dfv01$X1))
# vl2 <-  gsub("\\)","",dfv01$X2)
# vl2 <-  as.numeric(gsub(" ","",vl2))
# df_csvs05$vl1 <- vl1
# df_csvs05$vl2 <- vl2
value.vec <- df_csvs05$value
# define columns to keep
ctke    <- c("smplNm",
    "sample",
    "machine",
    # "vl1",
    # "vl2" ,
    "value",
    "spcNm"
    )
#
df_csvs05 <- df_csvs05[ctke]

nospc <- length(unique(df_csvs05$spcNm))
nosmp <- length(unique(df_csvs05$smplNm))
# make range of intervals to use for horizontal lines in the plot
y_values_hline <- seq(2.5, (2*(nospc)+0.5), 2)
# make range of intervals to use for vertical lines in the plot
x_values_hline <- seq(3.5, (3*(nosmp)+0.5), 3)

# <- lapply(df_csvs05, as.numeric)
# replace zeroes with NAs to avoid getting a point in the bubble chart
# df_csvs05$vl1[(df_csvs05$vl1==0)] <- NA
# df_csvs05$vl2[(df_csvs05$vl2==0)] <- NA

#use 'tidyr::pivot_longer' to get a long version of the data frame
# nvlcolm <- colnames(df_csvs05)[grepl("vl",colnames(df_csvs05))]
# df_csvs05.1 <- df_csvs05 %>% tidyr::pivot_longer(cols=c(nvlcolm),
#                                                  names_to = "vlc", 
#                                                  values_to="vln")
df_csvs05.1 <- df_csvs05
# paste together to get a column with both species name and sample name
df_csvs05.1$spA.smplN <- paste0(df_csvs05.1$spcNm,".",df_csvs05.1$smplNm)
# match to get colour for evaluation of eDNA level at sample 
# location per species
df_csvs05.1$cl.f.evl <- df_qPCRm02$cl.f.evl[match(df_csvs05.1$spA.smplN,df_qPCRm02$spA.smplN)]
# add white for NAs
df_csvs05.1$cl.f.evl[is.na(df_csvs05.1$cl.f.evl)] <- "white"
# copy column
df_csvs05.1$eval.c1 <- df_csvs05.1$cl.f.evl
# copy data frame
df_csvs05 <- df_csvs05.1

df_csvs05.1$vln <- df_csvs05.1$value
df_csvs05.1$vlc <- df_csvs05.1$machine

df_csvs05$vln <- df_csvs05$value
df_csvs05$vlc <- df_csvs05$machine

#
gplot4 <- ggplot(data=df_csvs05.1, aes(smplNm, sample, 
                  color = vlc, size = log10(vln))) + 
  theme_classic() +
  geom_point(aes(colour=machine)) +
  scale_size(range = c(0.1, 10), name="log10(molecules DNA/uL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # use section below to get rectangel backgrounds - however
  # this will require a second fill scale - and I was not able to 
  # get this second color scale working
  # geom_tile(aes(height=Inf,fill=as.factor(smplNm),y=0.005),
  #           position=position_nudge(x=-0.05), show.legend=F) +
  # scale_fill_manual(values=alpha(clfH,0.04)) +
  # # # 
  scale_color_manual(values=c("springgreen4", "#E69F00")) +
  geom_hline(yintercept = y_values_hline, colour = "gray64") +
  labs(color='machine')+
  scale_y_discrete(limits=rev) +
  xlab("MST sample") + ylab("species")
# see the ggplot object
#gplot4

# set parameter to check whether to save plot or not
bSaveFigures <- T
# make a filename with path
fnm02 <- paste0(wd00_01,"/Fig05_v02_compare_all_ddPCR_and_qPCR.png")
#p
if(bSaveFigures==T){
  ggsave(gplot4,file=fnm02,
         #width=297*0.6,height=210*1.2,
         height=297*0.6,width=210,
         units="mm",dpi=300)
}
#
df_csvs05$eDN.lv1 <- df_csvs05$vln
df_csvs05$eDN.lv1[is.na(df_csvs05$eDN.lv1)] <- 0
df_csvs05$eDN.lv1 <- as.numeric(df_csvs05$eDN.lv1)
# copy columns
df_csvs01$spcAbbr <- df_csvs01$speciesabbr
df_csvs01$f2 <- df_csvs01$sp.f2
df_csvs01$mcp_ddPCR_lod <- df_csvs01$ddlod
# get unique values for the species abbreviations
# for the LOD calculated per ddPCR analysis
df_qlod.ddP01 <- df_csvs01 %>% dplyr::group_by(spcAbbr) %>%
  distinct(spcAbbr, mcp_ddPCR_lod,f2,qlod)
# use match to transfer the qlod and factor for difference between qPCR
# and ddPCR
df_csvs05$f2 <- df_qlod.ddP01$f2[match(df_csvs05$spcNm,
                                       df_qlod.ddP01$spcAbbr)]
df_csvs05$dd.qlod <- df_qlod.ddP01$qlod[match(df_csvs05$spcNm,
                                              df_qlod.ddP01$spcAbbr)]
df_csvs05$mcp_ddPCR_lod <- df_qlod.ddP01$mcp_ddPCR_lod[match(
  df_csvs05$spcNm,df_qlod.ddP01$spcAbbr)]
# calculate loq10 levels of eDNA levels. #
# Add 1 to ensure log10 to zero is not infinite
df_csvs05$eDN.lv1.l10 <- log10(df_csvs05$eDN.lv1+1)
# now remove all lgo10 to 1, as these are zero, and are not needed on the plot
df_csvs05$eDN.lv1.l10[df_csvs05$eDN.lv1.l10==0] <- NA
# subset the dataframe to only comprise ddPCR findings
df_dP06 <- df_csvs05[(df_csvs05$machine=="ddPCR"),]
# copy the column with evaluation for colors
df_dP06$eval.c2 <- df_dP06$eval.c1

# evaluate a color coding for eDNA levels on ddPCR in relation to the qlod 
# and mcp_ddPCR_lod
df_dP06$eval.c2[df_dP06$eDN.lv1==0] <- "white"
df_dP06$eval.c2[df_dP06$eDN.lv1!=0 & df_dP06$eDN.lv1<df_dP06$dd.qlod] <- "yellow"
df_dP06$eval.c2[df_dP06$eDN.lv1>df_dP06$dd.qlod] <- "black"
# append back to combined data frame
df_csvs05$eval.c1[(df_csvs05$machine=="ddPCR")] <- df_dP06$eval.c2
# make a path for wrting a temporary file
wd00_01_tmpf05 <- paste0(wd00_01,"/tmp_df_csvs05.csv")
# write a temporary file for making a reprex
write.csv(df_csvs05, wd00_01_tmpf05)
# copy the evaluation column
df_csvs05$eval.c4 <- df_csvs05$eval.c1
#modify the evaluation category names
df_csvs05$eval.c1[df_csvs05$eval.c1=="white"] <- "no Cq" # no eDNA
df_csvs05$eval.c1[df_csvs05$eval.c1=="yellow"] <- "bLD" # below LOD
df_csvs05$eval.c1[df_csvs05$eval.c1=="orange"] <- "aLDbLQ" # above LOD below LOQ
df_csvs05$eval.c1[df_csvs05$eval.c1=="red"] <- "1aLQ" # 1 above LOQ
df_csvs05$eval.c1[df_csvs05$eval.c1=="black"] <- "aaLQ" # all above LOQ
# count the number of assays tested
no_of_as_Tst <- length(unique(df_csvs05$spcNm))
# ensure the plot library is loaded
library(ggplot2)
# begin the plot
plt_01a <- ggplot(df_csvs05, aes(smplNm, 
                                 sample, 
                                 size = eDN.lv1.l10,
                                 color = machine, 
                                 fill = eval.c1)) +
  theme_classic() +
  geom_point(shape = 21, stroke = 1.2) +
  geom_hline(yintercept = y_values_hline, colour = alpha("gray64",0.6)) +
  geom_vline(xintercept = x_values_hline, colour = alpha("gray64",0.6)) +
  scale_color_manual(values=alpha(c("springgreen4", "#E69F00"),0.9)) +
  scale_fill_manual(values=alpha(
    c("red","black","orange","yellow","white")
    ,0.5),name="evaluation") +
  scale_size(range = c(0.001, 10), 
             name="log10(molecules\n DNA/uL)") +
  scale_y_discrete(name="assay", 
                   breaks=arter_breaks, labels=arter_labels) +
  scale_x_discrete(name="watersample") +
  coord_flip() +
  labs(color='machine') +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1)) 
# adjust the size of the header categories on the legend key text
plt_01a  <- plt_01a + theme(legend.title=element_text(size=10))
#make a file to save to
fnm02 <- paste0(wd00_01,"/Fig05_v03_compare_all_ddPCR_and_qPCR.png")
#evaluate if the plot should be saved
if(bSaveFigures==T){
  ggsave(plt_01a,file=fnm02,
         #width=297*0.6,height=210*1.2,
         height=297*0.6,width=210,
         units="mm",dpi=300)
}

setwd(wd00_01)
# --------- regression ----------------------
# make a list of species abbreviations 
ls.spcAbbr <-  unique(df_csvs01$spcAbbr)
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
  # subset data frame based on species abbreviation
  sbs_csvs03 <- df_csvs01[df_csvs01$spcAbbr==spclat,]
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
write.csv(df_lm.mod02, file=paste0(wd00_01,"/table_w_linear_model_regression_results_v02.csv"))
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

write.csv(df_fit2, file=paste0(wd00_01,"/table_w_linear_model_regression_results_v03.csv"))
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
write.csv(df_fit2, file=paste0(wd00_01,"/table_w_linear_model_regression_results_v04.csv"))

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

#p2
ggsave(p2,filename="Fig06_v02_regr_ddPCR_qPCR.png",
       units="mm",width=210,height=297*0.8,dpi=300)

#subset to only comprise MST samples
dfMST1 <- df_csvs01[df_csvs01$smplTp=="MST",]
dfMST1 <- dfMST1 %>% dplyr::rename(ddPCR=mcp_ddPCR, qPCR=mcp_qPCR)
ls.spcAbbr2 <- unique(dfMST1$spcAbbr)
# iterate over species abbreviations to sum up in text the results
# start iteration
for (s in ls.spcAbbr2){
  # subset for species
  sb.dfMST1 <- dfMST1[dfMST1$spcAbbr==s,]
  
  sb.dfMST2 <- sb.dfMST1[sb.dfMST1$ddPCR>0,]
  MST.all.d <-  paste(sb.dfMST2$smplNm,collapse = ", ")
  n.of.all.ddPCR <- length(sb.dfMST2$smplNm)
  #print(paste0(s," - qPCR"))
  sb.dfMST3 <- sb.dfMST1[sb.dfMST1$qPCR>0,]
  n.of.all.qPCR <-length(sb.dfMST3$smplNm)
  MST.all.q <- paste(sb.dfMST3$smplNm,collapse = ", ")
  #only in qPCR
  onqP <- sb.dfMST3$smplNm[!sb.dfMST3$smplNm %in% sb.dfMST2$smplNm]
  #print(paste0(s," - only qPCR"))
  MST.only.q <- paste(onqP,collapse = ", ")
  n.of.only.q <- length(onqP)
  #only in qPCR
  ondP <- sb.dfMST2$smplNm[!sb.dfMST2$smplNm %in% sb.dfMST3$smplNm]
  #print(paste0(s," - only ddPCR"))
  MST.only.d <-  paste(ondP,collapse = ", ")
  n.of.only.d <- length(ondP)
  # sum up which are in both      
  #print(paste0(s," - in both"))
  inbth <- sb.dfMST3$smplNm[sb.dfMST3$smplNm %in% sb.dfMST2$smplNm]
  MST.inb <-  paste(inbth,collapse = ", ")
  n.of.inb <- length(inbth)
  #percentage in both
  prcinb <- round(n.of.inb/(n.of.only.d+n.of.only.q+n.of.inb)*100,0)
  txt <- paste0("For ",n.of.inb," prøver (",MST.inb,
                ") kunne både qPCR og ddPCR spore eDNA fra ",s," (Figur x)",
                ". Medens der med qPCR alene blev fundet eDNA i ",n.of.only.q,
                " vandprøver (",MST.only.q,
                "), og med ddPCR alene blev der fundet eDNA i ",n.of.only.d,
                " vandprøver (",MST.only.d,
                "). Dermed kunne begge platforme finde positivt eDNA signal i de samme ",prcinb,
                "% af alle positive sporinger ")
  print(txt)
}
# --------- alternative bubble (1) ----------------------

clfH <- rep(c("grey80","white"),length(unique(df_csvs05$smplNm)))

#Begin the plot
gplot4a <- ggplot(data=df_csvs05, 
                  aes(smplNm, sample, size = eDN.lv1.l10)) + 
  geom_tile(aes(height=Inf,fill=as.factor(smplNm),y=0.005),
            position=position_nudge(x=-0.05), show.legend=F) +
  scale_fill_manual(values=alpha(clfH,0.04)) +
  theme_classic() +
  geom_point(aes(colour=machine)) +
  scale_size(range = c(0.1, 10), name="log10(molecules DNA/uL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +

  scale_color_manual(values=c("springgreen4", "#E69F00")) +
  geom_hline(yintercept= y_values_hline, colour = "gray64") +
  labs(color='machine')+
  scale_y_discrete(limits=rev) +
  xlab("MST sample") + ylab("Species")
# see the ggplot object
#gplot4a

ggsave(gplot4a,filename="Fig07_v02_ddPCR_qPCR_1a.png",
       units="cm",width=16,height=8,dpi=300,scale=2)

# --------- alternative bubble (2) flip axes ----------------------

arter <-sort(unique(df_csvs05$sample))
arter_breaks <- arter[seq(1,length(arter),2)]
arter_labels <- stringr::str_remove(arter_breaks,"_ddPCR")
# begin a 2nd version of the plot
gplot4b <- ggplot(data=df_csvs05, aes(smplNm, sample, 
                                      size = eDN.lv1.l10)) + 
  geom_tile(aes(height=Inf,fill=as.factor(smplNm),y=0.005),
            position=position_nudge(x=-0.05), show.legend=F) +
  scale_fill_manual(values=alpha(clfH,0.04)) +
  theme_classic() +
  geom_point(aes(colour=machine)) +
  scale_size(range = c(0.1, 10), name="log10(molecules DNA/uL)") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust=1)) +

  scale_color_manual(values=c("springgreen4", "#E69F00")) +
  geom_hline(yintercept= y_values_hline, colour = "gray64") +
  labs(color='machine')+
  xlab("MST sample") + ylab("species") +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  scale_y_discrete(name="species", breaks=arter_breaks, 
                   labels=arter_labels)

#gplot4b

ggsave(gplot4b,filename="Fig08_v02_ddPCR_qPCR_1b.png",units="cm",
       width=12,height=12,dpi=300,scale=2)


# ------------------- bubble plot ---------------------------

df_csvs04 <- df_csvs03 %>%
  tidyr::pivot_longer(cols=c("mcp_qPCR","mcp_ddPCR"),
                      names_to="machine",
                      values_to="value") %>%
  mutate(value=ifelse(value==0,NA,value)) %>%
  mutate(machine=stringr::str_remove(machine,"mcp_"))
# make a vector for species assays to exclude because they
# did not amplify properly
not_ampl_prop <- c("Rhihar",
                   "Neomel",
                   "Oncgor")
# collapse the vector
not_ampl_prop <- paste(not_ampl_prop, collapse = "|" )
# use grepl to exlude the matches with the vector
df_csvs04.1 <- df_csvs04[!grepl(not_ampl_prop,df_csvs04$spcAbbr),]

df_csvs04 <- df_csvs04.1
# convert variables to factors
# machine
df_csvs04$machine <- factor(df_csvs04$machine,levels=c("qPCR","ddPCR"))
# sample
list_samples <- sort(unique(df_csvs04$smplNm))
df_csvs04$smplNm <- factor(df_csvs04$smplNm,levels=list_samples)
# species
list_species <- sort(unique(df_csvs04$spcAbbr))
df_csvs04$spcAbbr <- factor(df_csvs04$spcAbbr,levels=list_species)

# take every second sample
list_samples2 <- list_samples[seq(2,length(list_samples),2)]
#use this to assign a factor for shaded bars
df_csvs04 <- df_csvs04 %>%
  mutate(shade=ifelse(smplNm %in% list_samples2,1,0))
df_csvs04$shade <- factor(df_csvs04$shade,levels=c(0,1))
# add coordinates for vline to every species except the last one
df_csvs04 <- df_csvs04 %>%
  mutate(xline = ifelse(spcAbbr==list_species[length(list_species)],NA,2.5))
# transform eDNA levels to log10 scale , and add 1 for zero
df_csvs04$eDNA.lvl <- df_csvs04$value
df_csvs04$eDNA.lvl[is.na(df_csvs04$eDNA.lvl)] <- 0
df_csvs04$eDN.lv.l10 <- log10(df_csvs04$eDNA.lvl+1)
df_csvs04$eDN.lv.l10[(df_csvs04$eDN.lv.l10==0)] <- NA
#begin the plot
gplot4c <- ggplot(data=df_csvs04, aes(x=machine, y=smplNm)) + 
  geom_tile(aes(width=Inf,fill=shade), 
            alpha=0.3,show.legend=F) +
  theme_classic() +
  geom_point(aes(colour=machine, size=eDN.lv.l10),alpha=0.7) +
  scale_size(range = c(0.00001, 10), name="log10(molecules DNA/uL)") +
  
  theme(axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_fill_manual(values=c("grey64","white")) +
  scale_color_manual(values=c("#E69F00","springgreen4")) +
  #scale_color_manual(values=c("springgreen4","#E69F00")) +
  geom_vline(aes(xintercept=xline), colour = "gray64") +
  labs(color='machine')+
  scale_y_discrete(limits=rev) +
  ylab("MST sample") + xlab("Species") +
  # facet wrap allows for placing the abbreaviated assay name in the middle 
  # of each double column
  facet_wrap(.~spcAbbr, strip.position = "bottom",
             scales = "free_x",nrow=1) 
# see the ggplot object
#gplot4c

ggsave(gplot4c,filename="Fig08_v02_ddPCR_qPCR_1c.png",units="cm",
       width=12,height=12,dpi=300,scale=2)
#
size_labels <- seq(-2,6,2) 
size_breaks <- 10^(size_labels)
size_labels <- as.character(size_labels)
# get all unique species
unqspcAbbr <- unique(df_csvs04$spcAbbr)
nuspA <- length(unqspcAbbr)
# get the first 9 species
f8.1.spcAbr <- unqspcAbbr[1:9]
f8.2.spcAbr <- unqspcAbbr[10:nuspA]
# subset the list to only comprise the first 8 species
df_csvs04.01 <- df_csvs04[df_csvs04$spcAbbr  %in%  f8.1.spcAbr,]
df_csvs04.02 <- df_csvs04[df_csvs04$spcAbbr  %in%  f8.2.spcAbr,]
# make a plot with the first 9 species
gplot4d <- ggplot(data=df_csvs04.01, aes(x=machine, y=smplNm)) + 
  geom_tile(aes(width=Inf,fill=shade), 
            alpha=0.3,show.legend=F) +
  theme_classic() +
  geom_point(aes(colour=machine, size=eDN.lv.l10 ),alpha=0.7) +
  scale_size_area(trans="log1p",
                  name="log10(molecules DNA/uL)",
                  breaks=size_breaks,
                  labels=size_labels,
                  max_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_fill_manual(values=c("grey64","white")) +
  #scale_color_manual(values=c("springgreen4", "#E69F00")) +
  scale_color_manual(values=c("#E69F00","springgreen4")) +
  geom_vline(aes(xintercept=xline), colour = "gray64") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(color='machine')+
  scale_y_discrete(limits=rev) +
  ylab("MST sample") + xlab("Species") +
  facet_wrap(.~spcAbbr, strip.position = "bottom",
             scales = "free_x",nrow=1) 
# see plot
#gplot4d

ggsave(gplot4d,filename="Fig08_v02_ddPCR_qPCR_2a.png",units="cm",
       width=12,height=12,dpi=300,scale=2)
# make a plot for the last 9 species
gplot4d <- ggplot(data=df_csvs04.02, aes(x=machine, y=smplNm)) + 
  geom_tile(aes(width=Inf,fill=shade), 
            alpha=0.3,show.legend=F) +
  theme_classic() +
  geom_point(aes(colour=machine, size=eDN.lv.l10),alpha=0.7) +
  scale_size_area(trans="log1p",
                  name="log10(molecules DNA/uL)",
                  breaks=size_breaks,
                  labels=size_labels,
                  max_size = 10) +

  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_fill_manual(values=c("grey64","white")) +
  #scale_color_manual(values=c("springgreen4", "#E69F00")) +
  scale_color_manual(values=c("#E69F00","springgreen4")) +
  geom_vline(aes(xintercept=xline), colour = "gray64") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(color='machine')+
  scale_y_discrete(limits=rev) +
  ylab("MST sample") + xlab("Species") +
  facet_wrap(.~spcAbbr, strip.position = "bottom",
             scales = "free_x",nrow=1) 
# see the plot
#gplot4d
# save the plot
ggsave(gplot4d,filename="Fig08_v02_ddPCR_qPCR_2b.png",units="cm",
       width=12,height=12,dpi=300,scale=2)

#_______________________________________________________________________________
# start - map ddPCR and qPCR results on a map
#_______________________________________________________________________________
#define external working directory
#wd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2020/NOVANA_proever_2018_2019/"
#extwd00 <- "/home/hal9000/Documents/Documents/NIVA_Ansaettelse_2021/MONIS6"
extwd00 <- paste0(wd00,"/","MONIS6_2021_data")
extwd03 <- "output01a_2021_assembl_xls_MST_proevetagn_skema"
extwd01.2 <- "output01b_2021_merged_MST_csv_file"
extwd02 <- "output02_merged_txtfiles_from_mxpro_for_MONIS6"
# define external infile
extinf01 <- "sample_locations_MST_2021.csv"
# paste path and input file together
pth.extinf01<- paste0(extwd00,"/",extwd01.2,"/",extinf01)
# read sample locations
df_sl01 <- read_csv(pth.extinf01)
# 'breddegrad' is Danish for latitude
# 'længdegrad' is Danish for longitude
# use columns to make new with shorter names
df_sl01$dec_lat <- df_sl01$position_breddegrad_lok_pos
df_sl01$dec_lon <- df_sl01$position_lngdegrad_lok_pos
# copy column with unique sample number
df_sl01$unq.smpl.no  <- df_sl01$Vandprvenummer_unikt_nummer_U_Pr_Nr
# substitute in sample number
df_csvs05$smplNm2 <- gsub("-","",df_csvs05$smplNm)
# match to get lat and long
df_csvs05$dec_lon <- df_sl01$dec_lon[match(df_csvs05$smplNm2,df_sl01$unq.smpl.no)]
df_csvs05$dec_lat <- df_sl01$dec_lat[match(df_csvs05$smplNm2,df_sl01$unq.smpl.no)]
# match to get sampling date
df_csvs05$smpl.date <- df_sl01$dato_dato_inds[match(df_csvs05$smplNm2,
                                                    df_sl01$unq.smpl.no)]
# exclude wrong latitude
# https://sparkbyexamples.com/r-programming/r-subset-data-frame-by-column-value/
# Subset Rows by Checking values on Multiple Columns
df_csvs05.1 <- df_csvs05[(df_csvs05$dec_lat < 60 & df_csvs05$dec_lat > 50),]
df_csvs05.1 <- df_csvs05.1[(df_csvs05.1$dec_lon < 19 & df_csvs05.1$dec_lon > 1),]
# remove if machine is NA
df_csvs05.1 <- df_csvs05.1[!is.na(df_csvs05.1$machine),]
# remove if dec_lat and dec_lon is NA
df_csvs05.1 <- df_csvs05.1[!is.na(df_csvs05.1$dec_lon),]
df_csvs05.1 <- df_csvs05.1[!is.na(df_csvs05.1$dec_lat),]
# split string by delimiter to get only species abbrevation
df_csvs05.1$spcAbr <- sapply(strsplit(df_csvs05.1$sample,"_"), "[[", 1)
# make the sample date a character
df_csvs05.1$smpl.date <- as.character(df_csvs05.1$smpl.date)
# split string by delimiter to get day , month and year
df_csvs05.1$smpl.day <- sapply(strsplit(df_csvs05.1$smpl.date,"-"), "[[", 3)
df_csvs05.1$smpl.mnt <- sapply(strsplit(df_csvs05.1$smpl.date,"-"), "[[", 2)
df_csvs05.1$smpl.yea <- sapply(strsplit(df_csvs05.1$smpl.date,"-"), "[[", 1)
#make the day , month and year numeric
df_csvs05.1$smpl.day <- as.numeric(df_csvs05.1$smpl.day)
df_csvs05.1$smpl.mnt <- as.numeric(df_csvs05.1$smpl.mnt)
df_csvs05.1$smpl.yea <- as.numeric(df_csvs05.1$smpl.yea)
# remove unneeded column
df_csvs05.1$sample <- NULL 
# use tidyr to spread out the data frame to get columns for each machine
df_csvs05.2 <- df_csvs05.1 %>% tidyr::pivot_wider(
              names_from=machine,
              values_from = vln)
# rename columns using dplyr
df_csvs05.2 <- dplyr::rename(df_csvs05.2, qP.cc = qPCR)
df_csvs05.2 <- dplyr::rename(df_csvs05.2, dP.cc = ddPCR)
# make letters for each species
ulsp <- unique(df_csvs05.2$spcAbr)
nspo2<- length(ulsp)
letsp <- LETTERS[1:nspo2]
df_Lsp <- as.data.frame(cbind(letsp,ulsp))
df_Lsp$Ls3 <- paste0(letsp,")   ",ulsp)
df_csvs05.2$latspc3 <- df_Lsp$Ls3[match(df_csvs05.2$spcAbr,df_Lsp$ulsp )]
df_csvs05.2$llatspc3 <- df_Lsp$letsp[match(df_csvs05.2$spcAbr,df_Lsp$ulsp )]
# modify NAs to 0, then take log10 plus 1, and finally make all 0 NA
df_csvs05.2$dP.cc[is.na(df_csvs05.2$dP.cc)] <- 0
df_csvs05.2$qP.cc[is.na(df_csvs05.2$qP.cc)] <- 0
df_csvs05.2$l10.dP.cc <- log10(df_csvs05.2$dP.cc+1)
df_csvs05.2$l10.qP.cc <- log10(df_csvs05.2$qP.cc+1)
df_csvs05.2$l10.dP.cc[(df_csvs05.2$l10.dP.cc==0)] <- NA
df_csvs05.2$l10.qP.cc[(df_csvs05.2$l10.qP.cc==0)] <- NA
# Subset to get spring and fall samples -  in 2 data frame
# https://sparkbyexamples.com/r-programming/r-subset-data-frame-by-column-value/
# Subset Rows by Checking values on Multiple Columns
# get sampling from month 7 and below
df_csvs05.2s <- df_csvs05.2[(df_csvs05.2$smpl.mnt<=7),]
# get sampling from month 8 and above
df_csvs05.2f <- df_csvs05.2[(df_csvs05.2$smpl.mnt>=8),]

#_Get library to make ggplot maps
library(rnaturalearth)
# get worldmap
sf.df.world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
#begin ggplot for spring samples
p05 <- ggplot(data = sf.df.world) +
  geom_sf(color = "black", fill = "azure3") +
  theme_classic() +
  #add points for qPCR
  geom_point(data = df_csvs05.2s, 
             aes(x = dec_lon-0.16, y = dec_lat,
                 #size=sqrt((l10.qP.cc)/pi)),
                 size=l10.qP.cc,
                 fill=eval.c1),
             color=alpha(c("#E69F00"),c(0.7)),
             shape=21, stroke = 1.2) +
  #add points for ddPCR
  geom_point(data = df_csvs05.2s, 
             aes(x = dec_lon+0.16, y = dec_lat,
                 #size=sqrt((l10.dP.cc)/pi)),
                 size=l10.dP.cc,
                 fill=eval.c1),
             color=alpha(c("springgreen4"),c(0.7)),
             shape=21,stroke = 1.2) +
  
  #set the color of the points
  #scale_color_manual(values=c("springgreen4", "#E69F00")) +
  #Arrange in facets
  ggplot2::facet_wrap(. ~ llatspc3 + spcAbr,
                      ncol = 3,
                      labeller = label_bquote(cols = .(llatspc3) ~ .(" ") ~ italic(.(spcAbr))) ) +
  # 
  scale_fill_manual(values=alpha(c("red","black","orange","yellow","white"),0.5),
                    name="evaluation") +
  #change axis labels
  xlab("longitude") + ylab("latitude") +
  #xlab("længdegrad") + ylab("breddegrad") +
  # place legend on bottom
  theme(legend.position = "top") +
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(7, 17),
                    ylim = c(54.4, 59.0), 
                    expand = FALSE)
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05 + labs(size='log10(molecules DNA/uL)')
#see the plot
#p05t
# save the plot
ggsave(p05t,filename="Fig09_v02_map_spring_ddPCR_qPCR_1b.png",
       units="mm",
       width=210,height=297,dpi=300,scale=1.4)


#begin ggplot for fall samples
#begin ggplot for spring samples
p05 <- ggplot(data = sf.df.world) +
  geom_sf(color = "black", fill = "azure3") +
  theme_classic() +
  #add points for qPCR
  geom_point(data = df_csvs05.2f, 
             aes(x = dec_lon-0.16, y = dec_lat,
                 #size=sqrt((l10.qP.cc)/pi)),
                 size=l10.qP.cc,
                 fill=eval.c1),
             color=alpha(c("#E69F00"),c(0.7)),
             shape=21, stroke = 1.2) +
  #add points for ddPCR
  geom_point(data = df_csvs05.2f, 
             aes(x = dec_lon+0.16, y = dec_lat,
                 #size=sqrt((l10.dP.cc)/pi)),
                 size=l10.dP.cc,
                 fill=eval.c1),
             color=alpha(c("springgreen4"),c(0.7)),
             shape=21,stroke = 1.2) +
  
  #set the color of the points
  #scale_color_manual(values=c("springgreen4", "#E69F00")) +
  #Arrange in facets
  ggplot2::facet_wrap(. ~ llatspc3 + spcAbr,
                      ncol = 3,
                      labeller = label_bquote(cols = .(llatspc3) ~ .(" ") ~ italic(.(spcAbr))) ) +
  # 
  # 
  scale_fill_manual(values=alpha(c("red","black","orange","yellow","white"),0.5),
                    name="evaluation") +
  #change axis labels
  xlab("longitude") + ylab("latitude") +
  #xlab("længdegrad") + ylab("breddegrad") +
  # place legend on bottom
  theme(legend.position = "top") +
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(7, 17),
                    ylim = c(54.4, 59.0), 
                    expand = FALSE)
#change the header for the legend on the side, 
#this must be done for both 'fill', 'color' and 'shape', to avoid 
#getting separate legends
p05t <- p05 + labs(size='log10(molecules DNA/uL)')
#see the plot
#p05t

#getwd()
# save the plot
ggsave(p05t,filename="Fig09_v02_map_fall_ddPCR_qPCR_1b.png",
       units="mm",
       width=210,height=297,dpi=300,scale=1.4)

#_______________________________________________________________________________
# make letters for each species
ulsp <- unique(df_csvs05.1$spcAbr)
nspo2<- length(ulsp)
letsp <- LETTERS[1:nspo2]
df_Lsp <- as.data.frame(cbind(letsp,ulsp))
df_Lsp$Ls3 <- paste0(letsp,")   ",ulsp)
df_csvs05.1$latspc3 <- df_Lsp$Ls3[match(df_csvs05.1$spcAbr,df_Lsp$ulsp )]
df_csvs05.1$llatspc3 <- df_Lsp$letsp[match(df_csvs05.1$spcAbr,df_Lsp$ulsp )]
#make copies of columns to modify
df_csvs05.1$dec_lon.m <- df_csvs05.1$dec_lon
df_csvs05.1$dec_lat.m <- df_csvs05.1$dec_lat
# modify the longitude
df_csvs05.1$dec_lon.m[df_csvs05.1$machine=="ddPCR"] <- 
  (df_csvs05.1$dec_lon.m[df_csvs05.1$machine=="ddPCR"]-0.22)
df_csvs05.1$dec_lon.m[df_csvs05.1$machine=="qPCR"] <- 
  (df_csvs05.1$dec_lon.m[df_csvs05.1$machine=="qPCR"]+0.22)
#remove any rows with NAs for copy counts
#df_csvs05.1 <- df_csvs05.1[!is.na(df_csvs05.1$value),]
df_csvs05.1 <- df_csvs05.1[!is.na(df_csvs05.1$vln),]
# Subset to get spring and fall samples -  in 2 data frame
# https://sparkbyexamples.com/r-programming/r-subset-data-frame-by-column-value/
# Subset Rows by Checking values on Multiple Columns
# get sampling from month 7 and below
df_csvs05.1s <- df_csvs05.1[(df_csvs05.1$smpl.mnt<=7),]
# get sampling from month 8 and above
df_csvs05.1f <- df_csvs05.1[(df_csvs05.1$smpl.mnt>=8),]
#_Get library to make ggplot maps
library(rnaturalearth)
# get worldmap
sf.df.world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

#begin ggplot for spring samples
p06 <- ggplot(data = sf.df.world) +
  geom_sf(color = "black", fill = "azure3") +
  theme_classic() +
  #add points for ddPCR
  geom_point(data = df_csvs05.1s, 
             aes(x = dec_lon.m, y = dec_lat.m,
                 size=eDN.lv1.l10,
                 fill=eval.c1,
                 color=machine),
             shape=21, stroke=1.2) +
  #change size range of points
  # scale_size_continuous(range = c(
  #   (sqrt(1E-7)/pi),
  #   (sqrt(1E3)/pi) )) +
  scale_size_continuous(range = c(1E-1, 10)) +
  #set the color of the points
  scale_color_manual(values=alpha(c("springgreen4", "#E69F00"),c(0.9))) +
  scale_fill_manual(values=alpha(c("red","black","orange","yellow","white"),0.5),
                    name="evaluation") +
  #Arrange in facets
  ggplot2::facet_wrap(. ~ llatspc3 + spcAbr,
                      ncol = 3,
                      labeller = label_bquote(cols = .(llatspc3) ~ .(") ") ~ italic(.(spcAbr))) ) +
  #change axis labels
  xlab("longitude") + ylab("latitude") +
  #xlab("længdegrad") + ylab("breddegrad") +
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  labs(size='log10(molecules DNA/uL)') +
  labs(color='machine') +
  # place legend on bottom
  theme(legend.position = "top") +
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(7, 17),
                    ylim = c(54.4, 59.0), 
                    expand = FALSE)

#see the plot
p06t <- p06
#p06t
#getwd()
# save the plot
ggsave(p06t,filename="Fig09_v02_map_spring_ddPCR_qPCR_2b.png",
       units="mm",
       width=210,height=297,dpi=300,scale=1.4)

#unique(df_csvs05.1f$eval.c1)
#begin ggplot for fall samples
p06 <- ggplot(data = sf.df.world) +
  geom_sf(color = "black", fill = "azure3") +
  theme_classic() +
  #add points for ddPCR
  geom_point(data = df_csvs05.1f, 
             aes(x = dec_lon.m, y = dec_lat.m,
                 size=eDN.lv1.l10,
                 fill=eval.c1,
                 color=machine),
             shape=21,stroke=1.2) +
  #change size range of points
  # scale_size_continuous(range = c(
  #   (sqrt(1E-7)/pi),
  #   (sqrt(1E3)/pi) )) +
  scale_size_continuous(range = c(1E-1, 10)) +
  #set the color of the points
  scale_color_manual(values=alpha(c("springgreen4", "#E69F00"),c(0.9))) +
  scale_fill_manual(values=alpha(c("red","black","orange","yellow","white"),0.5),
                    name="evaluation") +
  #Arrange in facets
  ggplot2::facet_wrap(. ~ llatspc3 + spcAbr,
                      ncol = 3,
                      labeller = label_bquote(cols = .(llatspc3) ~ .(") ") ~ italic(.(spcAbr))) ) +
  #change axis labels
  xlab("longitude") + ylab("latitude") +
  #xlab("længdegrad") + ylab("breddegrad") +
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  labs(size='log10(molecules DNA/uL)') +
  labs(color='machine') +
  # place legend on bottom
  theme(legend.position = "top") +
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(7, 17),
                    ylim = c(54.4, 59.0), 
                    expand = FALSE)

#see the plot
p06t <- p06
#getwd()
# save the plot
ggsave(p06t,filename="Fig09_v02_map_fall_ddPCR_qPCR_2b.png",
       units="mm",
       width=210,height=297,dpi=300,scale=1.4)
#______________________________________________________________________________
# map in facet wrap by 1st and 2nd season 
# where 1 st season = jan-jul
# and 2nd season = aug-dec
#______________________________________________________________________________
# add an empty column
df_csvs05.1$season <- NA
# Subset Rows by Checking values on Multiple Columns
# get sampling from month 7 and below
df_csvs05.1$season[(df_csvs05.1$smpl.mnt<=7)] <- "Mar-Jun"
# get sampling from month 8 and above
df_csvs05.1$season[(df_csvs05.1$smpl.mnt>=8)] <- "Jul-Nov"

# get sampling from month 7 and below
df_csvs05.1$season[(df_csvs05.1$smpl.mnt<=7)] <- "1. Mar-Jun"
# get sampling from month 8 and above
df_csvs05.1$season[(df_csvs05.1$smpl.mnt>=8)] <- "2. Jul-Nov"


# paste together season and species abbreviation
df_csvs05.1$spcAbr.season <- paste0(df_csvs05.1$spcAbr," - ",df_csvs05.1$season)


# get worldmap
sf.df.world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
df_csvs05.1c <- df_csvs05.1[df_csvs05.1$spcAbr!="Corcas",]
unspA <- unique(df_csvs05.1c$spcAbr)
#unspA <- unspA[order(unspA)]
NunspA <- length(unspA)
NunspA.h <- floor(NunspA/2)
f8spcAb01 <- unspA[1:NunspA.h]
f8spcAb02 <- unspA[(NunspA.h+1):NunspA]
df_csvs05.1c01 <- df_csvs05.1c[(df_csvs05.1c$spcAbr  %in%   f8spcAb01) ,]
df_csvs05.1c02 <- df_csvs05.1c[(df_csvs05.1c$spcAbr  %in%   f8spcAb02),]
# reorder data frame by species name
df_csvs05.1c02 <- df_csvs05.1c02[order(df_csvs05.1c02$spcNm),]
#begin ggplot for samples arranged in facet wrap by season and species
p07 <- ggplot(data = sf.df.world) +
  geom_sf(color = "black", fill = "azure3") +
  theme_classic() +
  #add points for ddPCR
  geom_point(data = df_csvs05.1c01, 
             aes(x = dec_lon.m, y = dec_lat.m,
                 size=eDN.lv1.l10,
                 fill=eval.c1,
                 color=machine),
             shape=21,stroke=1.2) +
  #change size range of points
  # scale_size_continuous(range = c(
  #   (sqrt(1E-7)/pi),
  #   (sqrt(1E3)/pi) )) +
  scale_size_continuous(range = c(1E-1, 10)) +
  #set the color of the points
  scale_color_manual(values=alpha(c("springgreen4", "#E69F00"),c(0.9))) +
  scale_fill_manual(values=alpha(c("black","orange","yellow","white"),0.5),
                    name="evaluation") +
  #Arrange in facets
  ggplot2::facet_wrap(. ~ spcAbr.season,
                      ncol = 4,
                      labeller = label_bquote(cols = .(spcAbr.season)))  +
  #change axis labels
  xlab("longitude") + ylab("latitude") +
  #xlab("længdegrad") + ylab("breddegrad") +
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  labs(size='log10(molecules \nDNA/uL)') +
  labs(color='machine') +
  # place legend on bottom
  theme(legend.position = "top") +
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(7, 17),
                    ylim = c(54.4, 59.0), 
                    expand = FALSE)
# make legend across 2 rows
# see: https://stackoverflow.com/questions/27130610/legend-on-bottom-two-rows-wrapped-in-ggplot2-in-r
p07 <- p07 + guides(fill=guide_legend(nrow=2,byrow=TRUE),
                    size=guide_legend(nrow=2,byrow=TRUE),
                    color=guide_legend(nrow=2,byrow=TRUE))
#see the plot
p07t <- p07
#p07t

#getwd()
# save the plot
ggsave(p07t,filename="Fig10_v02_map_seasons_ddPCR_qPCR_3a.png",
       units="mm",
       width=210*0.8,height=297*0.6,dpi=300,scale=1.4)

#unique(df_csvs05.1c02$eval.c1)
#begin ggplot for samples arranged in facet wrap by season and species
p07 <- ggplot(data = sf.df.world) +
  geom_sf(color = "black", fill = "azure3") +
  theme_classic() +
  #add points for ddPCR
  geom_point(data = df_csvs05.1c02, 
             aes(x = dec_lon.m, y = dec_lat.m,
                 size=eDN.lv1.l10,
                 fill=eval.c1,
                 color=machine),
             shape=21, stroke=1.2) +
  #change size range of points
  # scale_size_continuous(range = c(
  #   (sqrt(1E-7)/pi),
  #   (sqrt(1E3)/pi) )) +
  scale_size_continuous(range = c(1E-1, 10)) +
  #set the color of the points
  scale_color_manual(values=alpha(c("springgreen4", "#E69F00"),c(0.9))) +
  scale_fill_manual(values=alpha(c("red","black","orange","yellow","white"),0.5),
                    name="evaluation") +
  
  #Arrange in facets
  ggplot2::facet_wrap(. ~ spcAbr.season,
                      ncol = 4,
                      labeller = label_bquote(cols = .(spcAbr.season)))  +
  #change axis labels
  xlab("longitude") + ylab("latitude") +
  #xlab("længdegrad") + ylab("breddegrad") +
  #change the header for the legend on the side, 
  #this must be done for both 'fill', 'color' and 'shape', to avoid 
  #getting separate legends
  #labs(size='log10(molekyler \nDNA/uL)') +
  labs(size='log10(molecules \nDNA/uL)') +
  labs(color='machine') +
  # place legend on bottom
  theme(legend.position = "top") +
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(7, 17),
                    ylim = c(54.4, 59.0), 
                    expand = FALSE)

# make legend across 2 rows
# see: https://stackoverflow.com/questions/27130610/legend-on-bottom-two-rows-wrapped-in-ggplot2-in-r
p07 <- p07 + guides(fill=guide_legend(nrow=2,byrow=TRUE),
                    size=guide_legend(nrow=2,byrow=TRUE),
                    color=guide_legend(nrow=2,byrow=TRUE))
#see the plot
p07t <- p07

# save the plot
ggsave(p07t,filename="Fig10_v02_map_seasons_ddPCR_qPCR_3b.png",
       units="mm",
       width=210*0.8,height=297*0.6,dpi=300,scale=1.4)
#_______________________________________________________________________________
# end - map ddPCR and qPCR results on a map
#_______________________________________________________________________________
