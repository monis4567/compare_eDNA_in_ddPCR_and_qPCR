#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
# comparing ddpcr and qpcr

# The code is based on the code provided by Gledis Guri
# Quantifying the detection sensitivity and precision of qPCR and ddPCR 
# mechanisms for eDNA samples
# Author
# 
# Gled Guri
# Published
# from the scientific paper:
# Guri, G., Ray.J.L., Shelton, A.O., Kelly, R.P., Præbel, K.,Allan, E.A., Yoccoz, N., Johansen, T., Wangensteen, O.S., Hanebrekke, T., Westgaard, J.-I., 2024. Quantifying the detection sensitivity and precision of qPCR and ddPCR mechanisms for eDNA samples. Ecology and Evolution, 14:e70678. https://doi.org/10.1002/ece3.70678
# 
# November 10, 2024
# https://html-preview.github.io/?url=https://github.com/gledguri/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples/blob/main/Code/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples.html
library(here)
# Load functions

extract_list_param <- function(stanMod){
  l <- stanMod@model_pars
  x <- list(extract_param(stanMod,"lp__")) %>% setNames("lp__")
  for (i in l) {
    x <- c(x,list(extract_param(stanMod,i)) %>% setNames(i))
  }
  x <- x[-which(names(x) == c("lp__"))]
  return(x)
}

cloglog <- function(theta) log(-log(1 - theta))

extract_param <- function(model=stanMod,parmeter="alpha"){
  return(summary(model, par = parmeter)$summary %>% 
           unlist()%>%as.data.frame%>%round(.,2))
}

inv.cloglog <- function(theta) 1-exp(-exp(theta))

logreg <- function(x,b0,b1,i) y=(1/(1+exp(-(b0+(b1*(x+i))))))

scientific_10 <- function(x) {
  c <- scales::scientific_format()(x)
  t <- gsub("1e", "10^", c)
  t2 <- gsub("10\\^\\+", "10\\^", t)
  str2expression(t2)}

# Load libraries

library(ggplot2)
library(dplyr)

#install.packages("rstudioapi")
#install.packages("rstan")
library(rstan)
#install.packages("here")
library(here);options(mc.cores = parallel::detectCores())
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ggplot2")
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(plyr) # Load plyr package
library(readr) # Load readr package
# Load data

getwd()

# make an object with the path to the data
pth_gh <-        "https://raw.githubusercontent.com/gledguri/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples/d8c6b248cf8b0c57127b9005cc47ec686099fc6c/Data/"
pf01 <- paste0(pth_gh,"/ddPCR_environmental_samples.csv")
# list all files on the github data webpage
f01 <- "Sample_metadata.csv"
f02 <- "ddPCR_environmental_samples.csv"
f03 <- "ddPCR_standard_samples.csv"
f04 <- "ddpcr_amplitude_cod_herring_standard_samples.csv"
f05 <- "ddpcr_amplitude_cod_saithe_standard_samples.csv"
f06 <- "qPCR_environmental_samples.csv"
f07 <- "qPCR_opt_Cod_Herring.csv"
f08 <- "qPCR_opt_Cod_Saithe.csv"
f09 <- "qPCR_standard_samples.csv"
f10 <- "sample_idx.csv"
# make a list of the files
lrf <-list(f01,f02,f03,f04,f05,f06,f07,f08,f09,f10)
# iterate over the list of files and read them into the environment
for (i in lrf){
  assign(i,read.csv(paste0(pth_gh,i)))
}
# Load data

# e_ddpcr <- read.csv(here('Data','ddPCR_environmental_samples.csv'))
# e_qpcr <- read.csv(here('Data','qPCR_environmental_samples.csv'))
# st_ddpcr <- read.csv(here('Data','ddPCR_standard_samples.csv'))
# st_qpcr <- read.csv(here('Data','qPCR_standard_samples.csv'))

e_ddpcr <- ddPCR_environmental_samples.csv
e_qpcr <- qPCR_environmental_samples.csv
st_ddpcr <- ddPCR_standard_samples.csv
st_qpcr <- qPCR_standard_samples.csv
#View(st_qpcr)
#View(st_ddpcr)

st_qpcr_cm <- st_qpcr %>% dplyr::filter(!is.na(Ct))
e_qpcr_cm <- e_qpcr %>% dplyr::filter(!is.na(Ct))

# sample_metadata <- read.csv(here('Data','Sample_metadata.csv'))
# stidx <- read.csv(here('Data','sample_idx.csv'))

sample_metadata <- Sample_metadata.csv
stidx <- sample_idx.csv

# View(stidx)
# View(sample_metadata)
# View(e_ddpcr)
# View(e_qpcr)

# Get the working directory path
wd00 <- getwd()
# read in qpcr data file from merged data csv files
dtdir <- "data/MONIS6_2021_data/output02_merged_txtfiles_from_mxpro_for_MONIS6"
flNm <- "outfile02_merged_mxpro_csvfls_MONIS6.csv"
indir <- paste0(wd00,"/",dtdir)
# read in the merged data file
qpcr_data01 <- read.csv(paste0(indir,
                               "/",
                               flNm),sep = ";",stringsAsFactors = F)

# make columns that has names that match what is used in the
# e_qpcr and e_ddpcr dataframes
qpcr_data01$Conc <- qpcr_data01$Quantitycopies
qpcr_data01$Species <- qpcr_data01$speciesabbr
qpcr_data01$Sample_Name <- qpcr_data01$smpltp
qpcr_data01$Ct <- qpcr_data01$CtdRn
# copy the data to a new dataframe
df_qpcr <- qpcr_data01

# read in the table with meatadata for the samples
indir03 <- paste0(wd00,"/data/MONIS6_2021_data/")
# make an object that holds the files in the directory
lst.fmtd <- list.files(indir03)
# get only the files that have the string "table04" in the name
fls_fmtd <- lst.fmtd[grepl("table04",lst.fmtd)]

# paste together the path to the file
fltbl04 <- paste0(indir03,"/",fls_fmtd)
# read in the file
df_metadata_tbl04 <- read.csv(fltbl04,sep = ";",stringsAsFactors = F)
library(dplyr)
# join data frames by the qpcrno and well columns
df_qpcr02 <- left_join(df_qpcr,df_metadata_tbl04,by=c("qpcrno",
                                                      "Well")) %>% 
  dplyr::select(-c(qpcrno,Well) ) 

# make a path for the directory that holds the ddpcr files
indir02 <- paste0(wd00,
                  "/data/data_ddpcr_runs")
# make a list of the files in the directory
lst.fl.id02 <- list.files(indir02)
# get only the files that have the string "dddpcr" in the name
fls_ddpcr <- lst.fl.id02[grepl("ddpcr",lst.fl.id02)]
fls_ddpcr <- fls_ddpcr[grepl("csv",fls_ddpcr)]
# exclude files that have "0020" in the name
fls_ddpcr <- fls_ddpcr[!grepl("0020",fls_ddpcr)]
fls_ddpcr <- fls_ddpcr[!grepl("0030",fls_ddpcr)]
# add the path to the file names
lst_fpthinf <- paste0(indir02,"/",fls_ddpcr)
# Read in all csv files with ddpcr results
names(lst_fpthinf) <- basename(fls_ddpcr)
# To read in all files
combined_data <- lst_fpthinf %>%
  lapply(read_csv, col_types=cols(.default = col_character())) %>%
  bind_rows(.id="file")
# rename the resulting dataframe
df_csvs01 <- combined_data

# split the string by delimiter
dts <- data.frame(do.call('rbind',
                          strsplit(as.character(df_csvs01$file)
                                   ,'_',fixed=TRUE)))
#View(dts)
# assign the ddpcr no and species abbreviation to the dataframe
df_csvs01$ddpcrno <- dts$X1
df_csvs01$scpabbr <- dts$X2
# replace the string "Myaara" with "Myaare" in the scpabbr column
df_csvs01$scpabbr <- gsub("Myaara","Myaare",df_csvs01$scpabbr)
unique(df_csvs01$scpabbr)
unique(df_csvs01$ddpcrno)

# read in multiple xls files with the setup for ddpcr
lst_xls_ddpcr <- list.files(paste0(wd00,"/data/data_ddpcr_runs"),
                            pattern = "*.xls",full.names = T)
# subset in the list of xls-files that have "setup" in the name
lst_xls_ddpcr <- lst_xls_ddpcr[grepl("plate_setup",lst_xls_ddpcr)]
# and exclude the files that have 'stddiltest' in the name
lst_xls_ddpcr <- lst_xls_ddpcr[!grepl("stddiltest",lst_xls_ddpcr)]
# substitue the path to the file names, to get only the file name
lst_flnm_xls <- sub(".*/", "", lst_xls_ddpcr)

library(readxl)
# Read in all csv files with ddpcr results
names(lst_xls_ddpcr) <- basename(lst_flnm_xls)
# To read in all files
combined_data <- lst_xls_ddpcr %>%
# se this example to read in multiple files
# https://stackoverflow.com/questions/13441204/using-lapply-and-read-csv-on-multiple-files-in-r
    lapply(function(i){
    readxl::read_excel(i,skip=25, col_types = c("text")) 
  }) %>%
  bind_rows(.id="file")
# rename the resulting dataframe
df_xls01 <- combined_data
# split the string by delimiter
dts <- data.frame(do.call('rbind',
                          strsplit(as.character(df_xls01$file)
                                   ,'_',fixed=TRUE)))
#View(dts)
# assign the ddpcr no and species abbreviation to the dataframe
df_xls01$ddpcrno <- dts$X3
df_xls01$scpabbr <- gsub("(*)\\.xls","\\1",dts$X4)
df_xls01$Wellsample <- df_xls01$WellName
df_csvs01$WellNumber <- df_csvs01$Well
library(dplyr)
# combine the two dataframes by the ddpcrno and scpabbr columns
df_ddpcr <- dplyr::left_join(df_csvs01,df_xls01,by=c("ddpcrno",
                                                     "scpabbr",
                                                     "WellNumber")) %>% 
  dplyr::select(-file.x,-file.y)
# make common columns
df_ddpcr$MST_sample_number <- df_ddpcr$Wellsample
df_metadata_tbl04$MST_sample_number <- df_metadata_tbl04$MST.nummer
# subsitute in the column names
clnmm04 <- colnames(df_metadata_tbl04)
clnmm04 <- stringi::stri_trans_general(clnmm04, "latin-ascii")
clnmm04[grepl("ubit.",clnmm04)] <- "conc_qubit_ng_uL"
colnames(df_metadata_tbl04) <- clnmm04

# define the columns to keep
ckeep <- c(
"conc_qubit_ng_uL", "Dato.for.ekstraktion", 
           "Ekstraktionsnummer", "Noter2", "Elueringvolumen.AE.buffer.uL", 
           "Sub_Nmb_for_extraction", "Dato_for_indsml", "X.Noter2", "U_Pr_Nr", 
           "MSTno_di", "wdt", "Navn_Inds", "skib", "Inst", "Dato_inds", 
           "Lok_omr01", "lokalitet_vanda", "lok_pos_lat", "lok_pos_lon", 
           "max_Dyb", "Dyb_Str_m", "Vwf_mL", "vandsoejllagd", "Dyb_lagd01", 
           "dato_80C", "Temp_Inds", "Bemaerkn", "MST.n_Dinds", "MST.nummer", 
           "MST_STEXno", "smplNoSTEX", "MST_sample_number")
# keep only the columns that are in the ckeep vector
df_metadat05 <- df_metadata_tbl04[ckeep]
# only keep unique rows - see: https://stackoverflow.com/questions/68644516/how-to-keep-unique-rows-in-a-data-frame-using-dplyrdistinct-while-specifyi
df_metadat05 <- df_metadat05 %>% dplyr::distinct()
# omit the rows that have NA in the MST_sample_number column
df_metadat05 <- df_metadat05[!is.na(df_metadat05$MST_sample_number),]


# get the columns names, and find the column that
# has the string "description 1" in the name
# and substitute this column header with "smpl_descript_1"
clnm_ddpcr <- colnames(df_ddpcr)
idxsdescr1 <- which(grepl("description 1",clnm_ddpcr))
idxsdescr2 <- which(grepl("description 2",clnm_ddpcr))
idxsdescr3 <- which(grepl("description 3",clnm_ddpcr))
idxsdescr4 <- which(grepl("description 4",clnm_ddpcr))
# substitute the column names
clnm_ddpcr[idxsdescr1] <- "smpl_descript_1"
clnm_ddpcr[idxsdescr2] <- "smpl_descript_2"
clnm_ddpcr[idxsdescr3] <- "smpl_descript_3"
clnm_ddpcr[idxsdescr4] <- "smpl_descript_4"
colnames(df_ddpcr) <- clnm_ddpcr

# in the 'df_ddpcr$MST_sample_number' column, the 'std' samples
# needs to be with small letters, for the 'std' samples
# use gsub to replace the 'STD' with 'std' in 
# the 'MST_sample_number' column
df_ddpcr$MST_sample_number <- gsub("STD",
                                     "std",
                                     df_ddpcr$MST_sample_number)
# for the rows that have NAs for the MST_sample_number column
# then replace this with the value in the smpl_descript_1 column
df_ddpcr$MST_sample_number[is.na(df_ddpcr$MST_sample_number)] <- df_ddpcr$smpl_descript_1[is.na(df_ddpcr$MST_sample_number)]
df_ddpcr$MST_sample_number[(df_ddpcr$MST_sample_number=="0" & df_ddpcr$smpl_descript_1=="STD")] <- "std"
# check the unique values in the MST_sample_number column
# replace the '-' in the MST_sample_number column with an empty string
df_qpcr02$MST_sample_number <- gsub("-","",df_qpcr02$Sample_Name )
df_qpcr$MST_sample_number <- gsub("-","",df_qpcr$Sample_Name )
# copy the column to a new column, as the next line will overwrite the 
# original column
df_ddpcr$MST_sample_number2 <- df_ddpcr$MST_sample_number
# substitute the string "MST2021-" with "MST2021" in the MST_sample_number column
df_ddpcr$MST_sample_number <- gsub("MST2021-","MST2021",df_ddpcr$MST_sample_number)
# identify the samples with 3 variable
smplw3v <- df_ddpcr$MST_sample_number[grepl(".*_.*_.*",df_ddpcr$MST_sample_number)]
# split the string by the delimiter "_" 
# to get the standard dilution step
std.dil_step <- sapply(strsplit(smplw3v,"_"), "[[", 3)
std.dil_Nm <- sapply(strsplit(smplw3v,"_"), "[[", 2)
std.dil <- paste0(std.dil_Nm,"_",std.dil_step)
# add back the standard dilution to the MST_sample_number column
# to replace the values that have the species abbreviation included
df_ddpcr$MST_sample_number[grepl(".*_.*_.*",df_ddpcr$MST_sample_number)] <- std.dil
# get the rows that have the string "STD" in the MST_sample_number column
wllNm.STD <- df_ddpcr[grepl("std",df_ddpcr$MST_sample_number),]

# get the standard dilution levels from the Hemsan assay
stdillvl <- df_ddpcr$MST_sample_number[grepl("01",df_ddpcr$Well) & 
                  grepl("std",df_ddpcr$MST_sample_number) &
                  grepl("Hemsan",df_ddpcr$Target)]
# get the corresponding well names
WellNm <- df_ddpcr$Well[grepl("01",df_ddpcr$Well) & 
                  grepl("std",df_ddpcr$MST_sample_number) &
                  grepl("Hemsan",df_ddpcr$Target)]
# combine these columns into a dataframe
df_wllstd_Hemsan <- as.data.frame(cbind(stdillvl,WellNm))

# use the 'df_wllstd_Hemsan' to match the standard dilution levels
# for the samples that have the string "STD" in the MST_sample_number column
# and are missing the correct standard dilution level
wllNm.STD$smpl_descript_1 <- 
  df_wllstd_Hemsan$stdillvl[match(wllNm.STD$Well,
                                  df_wllstd_Hemsan$WellNm)]
# replace the values that only have STD
df_ddpcr[grepl("std",df_ddpcr$MST_sample_number),] <- wllNm.STD
# make an initial concentraion column
df_ddpcr$init_concentration <- NA
# use gsub on the 'smpl_descript_1' column 
# to get the initial concentration from only the 'std' samples
df_ddpcr$init_concentration[grepl("std",df_ddpcr$MST_sample_number)] <- 
  gsub(".*_","",
       df_ddpcr$smpl_descript_1[grepl("std",df_ddpcr$MST_sample_number)])
# make the 'init_concentration' column numeric
df_ddpcr$init_concentration <- as.numeric(df_ddpcr$init_concentration)

# get all the column names and place in a vector
clnm_dd <- colnames(df_ddpcr)
# get the index number for the column that has the string "Conc" in the name
idx_conc_clm <- which(grepl("Conc",clnm_dd))
# change the column name to "conc_copies_per_uL"
clnm_dd[idx_conc_clm] <- "conc_copies_per_uL"
# get the index number for the column that  has the string "Accepted Droplets" in the name
idx_conc_clm <- which(grepl("Accepted Droplets",clnm_dd))
# change the column name  to "Accepted_Droplets"
clnm_dd[idx_conc_clm] <- "Accepted_Droplets"
# replace the column names in the dataframe
colnames(df_ddpcr) <- clnm_dd
# the StanmodM2 needs a column that is called 'Tot_drop'
# which must be the same as the 'Accepted_Droplets' column
df_ddpcr$Tot_drop <- df_ddpcr$Accepted_Droplets

# select only the columns that are used in the 'e_ddpcr' dataframe
# that originally was provided in the code developed by Gledis Guri
# to end up with a data frame that resembles the input that is used in
# the next 'stan model part' in the code
df_ddpcr02 <- df_ddpcr %>% dplyr::select(Well, 
                    MST_sample_number,
                    scpabbr, 
                    conc_copies_per_uL,
                    Accepted_Droplets, 
                    Positives, 
                    Negatives,
                    Tot_drop,
                    init_concentration,
                    ddpcrno)
# only keep some of the columns
# is does not appear to matter if the column included are exact
# matches to the columns used in the example code from Gleidis Guri
df_metadat06 <- df_metadat05 %>% dplyr::select(Bemaerkn,
                        conc_qubit_ng_uL,
                        Dato.for.ekstraktion, 
                        dato_80C,
                        Dato_for_indsml,
                        Dato_inds, 
                        Dyb_lagd01, 
                        Dyb_Str_m, 
                        Ekstraktionsnummer, 
                        Elueringvolumen.AE.buffer.uL, 
                        Inst,
                        lokalitet_vanda, 
                        Lok_omr01, 
                        lok_pos_lat,
                        lok_pos_lon, 
                        max_Dyb, 
                        MST.nummer, 
                        MST.n_Dinds,
                        MSTno_di, 
                        MST_sample_number,
                    MST_STEXno,
                    Navn_Inds,
                    Noter2,
                    skib,
                    smplNoSTEX,
                    Sub_Nmb_for_extraction, 
                    Temp_Inds, 
                    U_Pr_Nr, 
                    vandsoejllagd,
                    Vwf_mL, 
                    wdt)
# get the sampling day, month and year
df_metadat06$smpl_day <- format(as.Date(df_metadat06$Dato_inds), "%d")
df_metadat06$smpl_mnt <- format(as.Date(df_metadat06$Dato_inds), "%m")
df_metadat06$smpl_yea <- format(as.Date(df_metadat06$Dato_inds), "%Y")
# convert these columns to numeric format
library(dplyr)
df_metadat06 <- df_metadat06 %>%
  dplyr::mutate(across(c(smpl_day,
                  smpl_mnt,
                  smpl_yea), as.numeric))
# subset the metadata to only comprise samples collected in 2021
# as some of the MST numbers unfortunately are not unique - as they
# were supposed to be
# Limiting the meta data to only comprise samples from 2021
# will help work around that the left join part following in the next section
# will not complain about the 'many-to-many relationship' 
# between the two dataframes
df_metadat06 <- df_metadat06[(df_metadat06$smpl_yea==2021),]
unique(df_metadat06$smpl_yea)
nrow(df_metadat06)

# combine ddpcr and metadata table by wellnumber and 
# MST sample number usin dplyr::left_join
df_ddpcr03 <- dplyr::left_join(df_ddpcr02,
                               df_metadat06,
                               by=c("MST_sample_number"))
  #                               %>% 
  # select(-c(MST_sample_number))


# combine ddpcr and metadata table by wellnumber and 
# MST sample number usin dplyr::left_join
df_qpcr03 <- dplyr::left_join(df_qpcr,
                               df_metadat06,
                               by=c("MST_sample_number")) 
  #                                %>% 
  # select(-c(MST_sample_number))
# replace the column name that is spelled wrong
clnm <-colnames(df_ddpcr03)
clnm[grepl("sc.abb",clnm)] <- "spcabbr"
colnames(df_ddpcr03) <- clnm
# replace the column name in the qpcr data frame that 
# has a different spelling for the "spcabbr"
clnm <-colnames(df_qpcr03)
clnm[grepl("sp.*abb",clnm)] <- "spcabbr"
colnames(df_qpcr03) <- clnm

spabb.dd <- unique(df_ddpcr03$spcabbr)
spabb.q <- unique(df_qpcr03$spcabbr)
# subset the qpcr data frame to only comprise the rows where the
# species abbreviations match what is occuring in the ddpcr
# data frame
# make a vector that can be used for subsetting the 'df_qpcr03' data frame
spc_abb_in.dd <- paste0(spabb.dd,collapse = "|")
# then subset the qpcr data frame so that the qpcr data frame
# only contains the species abbreviations that are included in 
# the ddpcr data frame
df_qpcr03 <- df_qpcr03[grepl(spc_abb_in.dd,df_qpcr03$spcabbr),]





# check out the contents of the st_ddpcr data frame
# that comes with the original code, to see what numbers are associated
# for each 'st_idx'.
st_ddpcr %>% dplyr::distinct(Species,Sample_name,st_idx) %>%
  dplyr::arrange(Species, Sample_name,st_idx) 
st_ddpcr %>% dplyr::distinct(Sample_name,st_idx) %>%
  dplyr::arrange(Sample_name,st_idx) 

st_ddpcr %>% dplyr::distinct(Sample_name,st_idx,int_concentation) %>%
  dplyr::arrange(st_idx,Sample_name,int_concentation) 

st_qpcr %>% dplyr::distinct(Sample_name,st_idx,int_concentation) %>%
  dplyr::arrange(st_idx,Sample_name,int_concentation) 

st_qpcr %>% dplyr::distinct(Species,Sample_name,st_idx) %>%
  dplyr::arrange(Species, Sample_name,st_idx) 
st_qpcr %>% dplyr::distinct(Sample_name,st_idx) %>%
  dplyr::arrange(Sample_name,st_idx) 
# # Gledis Guri explained in an email:
# # Regarding the column you mentioned
# pres = presence / absence of the qPCR observation. 
# It takes values of either 0 = absence or 1 = presence of qPCR amplification.
# species_idx = species index or assay index. 
# In my case I have 1, 2, and 3 for Cod, Herring, and Saithe respectively. 
# Since you have 16 assays you will be numbering them from 1 - 16.
# st_idx = station index or sample index. 
# So if you have 2 samples that you run 3 technical replicates in 
# qPCR for each sample then the indexing will be 1,1,1,2,2,2. 
# This indicates which row belongs to which sample. 
# If you have 10 samples and 3 tech replicates you will have a vector of 
# length 30 with numbers from 1 - 10.
# sp_st_idx = species - station index. This is a combination of the two 
# previous indexes (species_idx and st_idx). Say you have the 
# same 2 samples (each with 3 technical replicates) for each of the 3 assays. 
# The index vectors will be as follow:
#   species_idx [1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3]
# st_idx         [1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2]
# sp_st_idx    [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6] 
# Your qPCR observation will be 2 (samples) * 3 (tech rep) * 3 (assays) = 18 rows 
# therefore sp_st_idx will be a vector of length 18 with numbers from 1 - 6.

# for the 'e_qpcr' and the 'e_ddpcr' data frames
# the 'st_idx' column represents a unique forthrunning number
# that is associated with each sample
e_qpcr %>% dplyr::distinct(Species,Sample_name,st_idx) %>%
  dplyr::arrange(Species, Sample_name,st_idx) 
e_qpcr %>% dplyr::distinct(Species,species_idx) %>%
  dplyr::arrange(Species,species_idx) 
e_qpcr %>% dplyr::distinct(Sample_name,st_idx) %>%
  dplyr::arrange(Sample_name,st_idx) 
e_ddpcr %>% dplyr::distinct(Species,Sample_name,st_idx) %>%
  dplyr::arrange(Species, Sample_name,st_idx) 
e_ddpcr %>% dplyr::distinct(Sample_name,st_idx) %>%
  dplyr::arrange(Sample_name,st_idx)

# well number in the 'df_qpcr03' data frame does not have a zero in front of the number that 
# represents the column number of the plate. Use the 'str_pad' function to pad with
# zeroes, but first use 'gsub' to differentiate between the letters and the numbers in 
# the 'Well' column
# use gsub to differentiate between the letters and the numbers in the 'Well' column
Well.let <- gsub("[[:digit:]]","",df_qpcr03$Well)
Well.nmb <- gsub("[[:alpha:]]","",df_qpcr03$Well)
#pad with zeros to two characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
Well.nmb <- stringr::str_pad(Well.nmb, 2, pad = "0")
#paste together the two columns to get a Well column with preceeding zeroes in front of
# the column number of the plate
df_qpcr03$Well <- paste0(Well.let,Well.nmb)
# reorder the data frame
df_ddpcr03 <- df_ddpcr03 %>% dplyr::arrange(ddpcrno, Well, MST_sample_number, spcabbr) 
df_qpcr03 <- df_qpcr03 %>% dplyr::arrange(qpcrno, Well, MST_sample_number, spcabbr) 
# the concentration column in the 'df_qpcr03'  represents the Quantity copies
# this needs to be a numeric value, and the 'NoCt' values need to be replaced with 0
df_qpcr03$Conc[(df_qpcr03$Conc=="NoCt")] <- 0
df_qpcr03$Ct[(df_qpcr03$Ct=="NoCt")] <- 0
# then convert the 'Conc' and 'Ct' columns and the
# the smpl_day, smpl_mnt, and smpl_yea to numeric values
df_qpcr03 <- df_qpcr03 %>% dplyr::mutate_at(c('Conc',
                                              'Ct',
                                            'smpl_day',
                                            'smpl_mnt',
                                            'smpl_yea'), as.numeric)
# use the 'str' function to check out which are characters and which are numeric
df_qpcr03 %>% str()

# make a column with initial concentration, start by adding NAs, as the
# content of this column will depend on the concentration 
# level of the standard dilution series
df_qpcr03$init_concentration <- NA
# get all the rows that represent the standard dilution series
stds <- df_qpcr03$smpltp[grepl("std",df_qpcr03$smpltp)]
# get the standard dilution level, by using gsub to only get the numeric
# value of the dilution level
std.dil_lvl <- gsub("std","",stds)
# evaluate whether the sample is an NTC sample, in such a case
# the standard dilution level is 0, and should be modified accordingly
# check the values in the 'std.dil_lvl' column, using grepl, and grep
# for 'NTC'
std.dil_lvl[grepl("NTC",std.dil_lvl)] <- 0
# convert the 'std.dil_lvl' column to numeric
std.dil_lvl <- as.numeric(std.dil_lvl)
# use the same subsetting as above, grep'ing for 'std' in
# the 'smpltp' column, to get the rows that represent the standard dilution series
# the replace the 'init_concentration' column with the values in the 'std.dil_lvl' column 
df_qpcr03$init_concentration[grepl("std",df_qpcr03$smpltp)] <- std.dil_lvl

# the 'df_ddpcr03' data frame has the 
# 'Accepted_Droplets' column that is a character
# convert this to a numeric value
# and convert all other numeric columns
df_ddpcr03 <- df_ddpcr03 %>% dplyr::mutate_at(c('Accepted_Droplets',
                                                'Positives',
                                                'Negatives',
                                                'Tot_drop',
                                                'conc_qubit_ng_uL',
                                                'conc_copies_per_uL',
                                                'init_concentration',
                                                'smpl_day',
                                                'smpl_mnt',
                                                'smpl_yea'), as.numeric)
# use the 'str' function to check out which are characters and which are numeric
df_ddpcr03 %>% str()
# copy the 'Accepted_Droplets' column to a new column that is called 'Total_drop'
df_ddpcr03$Total_drop <- df_ddpcr03$Accepted_Droplets
# use gsub to replace the 'Myaara' with 'Myaare' in the 'spcabbr' column
df_ddpcr03$spcabbr <- gsub("Myaara","Myaare",df_ddpcr03$spcabbr)
# the df_ddpcr03 and the df_qpcr03 data frames need a column named
# 'Species', based on the species abbreviation
df_ddpcr03$Species <- df_ddpcr03$spcabbr
df_qpcr03$Species <- df_qpcr03$spcabbr
unique(df_ddpcr03$Species) %in% unique(df_qpcr03$Species)
# subset to only comprise the 'Hemtak' species
df_ddpcr03.Hemtak <- df_ddpcr03[grepl("Hemtak",df_ddpcr03$spcabbr),]
df_ddpcr03.NA_for_init_conc <- df_ddpcr03[is.na(df_ddpcr03$init_concentration),]
# get only the weel that are from the first column in the plate
df_ddpcr03.NA_for_init_conc <- df_ddpcr03.NA_for_init_conc[grepl("01",df_ddpcr03.NA_for_init_conc$Well),]


pth_to_gh_stc <- "https://github.com/gledguri/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples/tree/d8c6b248cf8b0c57127b9005cc47ec686099fc6c/Code"
#                  "https://github.com/gledguri/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples/blob/d8c6b248cf8b0c57127b9005cc47ec686099fc6c/Code/Stan_qPCR_ddPCR_prob_of_det.stan"
pth_to_gh_stc <- "https://raw.githubusercontent.com/gledguri/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples/refs/heads/main/Code/"
wd00 <- getwd()
#Create a directory to put the code obtained from github into
wd00_wdc <- paste0(wd00,"/code")
#Delete any previous versions of the output directory
unlink(wd00_wdc, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wdc)
#set the output directory
outdir01 <- wd00_wdc

#no_sbs <- "01"
#Create a directory to put the resulting figures into
wd00_wd02 <- paste0(wd00,
                    paste0("/output13_",#,no_sbs,
                           "_my_df_modified_to_match",
                           "_G_Guri_setup"
                    ))
#Delete any previous versions of the output directory
unlink(wd00_wd02, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wd02)
#set the output directory
outdir02 <- wd00_wd02

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# write out the 'df_ddpcr03' and the 'df_qpcr03' data frames
# as tsv files  
write.table(df_ddpcr03,
            file = paste0(wd00_wd02,"/df_ddpcr03.tsv"),
            sep = "\t",
            row.names = F)
write.table(df_qpcr03,
            file = paste0(wd00_wd02,"/df_qpcr03.tsv"),
            sep = "\t",
            row.names = F)

#
#read in the tsv tables
df_ddpcr03 <- read.table(paste0(wd00_wd02,"/df_ddpcr03.tsv"),
           sep = "\t",
           header = T) 
df_qpcr03 <- read.table(paste0(wd00_wd02,"/df_qpcr03.tsv"),
           sep = "\t",
           header = T)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# 
# 
# # make a column on the ddPCR data frame that can be used
# # for distinguishing the 'NTC', the 'Standard' and the 'eDNA' samples
# # make an empty column with NAs for the 'WellType' column
# df_ddpcr03$WellType <- NA
# # if the column 'MST_sample_number' has the string 'std' in it
# # the make it a 'Standard' sample
# df_ddpcr03$WellType[grepl("std",df_ddpcr03$MST_sample_number)] <- "Standard"
# # if the column 'MST_sample_number' has the string 'NTC' in it
# # the make it a 'NTC' sample
# df_ddpcr03$WellType[grepl("NTC",df_ddpcr03$MST_sample_number)] <- "NTC"
# # the rest of the samples are 'Unknown' samples
# df_ddpcr03$WellType[is.na(df_ddpcr03$WellType)] <- "Unknown"
# 
# # copy the 'conc_copies_per_uL' column to a 
# # new column that is called 'Conc'
# df_ddpcr03$Conc <- df_ddpcr03$conc_copies_per_uL
# #check if any has NAs in 'Conc' column and , if they have
# # then replace with 0
# df_ddpcr03$Conc[is.na(df_ddpcr03$Conc)] <- 0
# # copy the column named 'init_concentration' to a new column
# df_qpcr03$int_concentation <- df_qpcr03$init_concentration
# 
# # get the unique species abbreviations in the ddpcr data frame
# # and reorder them alphabetically
# spcabbr.d <- unique(df_ddpcr03$spcabbr)
# spcabbr.d <- spcabbr.d[order(spcabbr.d)]
# length(spcabbr.d)
# 
# 
# # get the unique species abbreviations in the qpcr data frame
# # and reorder them alphabetically
# spcabbr.q <- unique(df_qpcr03$spcabbr)
# spcabbr.q <- spcabbr.q[order(spcabbr.q)]
# # check that they have the same number of species
# length(spcabbr.d)
# length(spcabbr.q)
# # also check it is the same species
# # that are included in both data frames
# spcabbr.d %in% spcabbr.q
# spcabbr.q %in% spcabbr.d
# 
# #sbs_spcs <- spcabbr.q[1:2]
# spabb.q <- spabb.q[order(spabb.q)]
# #idxnsbs <- which(spabb.q %in% sbs_spcs)
# 
# # In the ggplots , later on, a color is required for each species
# # start by making a data frame
# # get one range of colors
# clr01 <- palette.colors(palette = "Okabe-Ito")
# # get another range of colors
# clr02 <- palette.colors(palette = "Polychrome")
# # combine the two ranges of colors
# clr03 <- c(clr01,clr02)
# # use the combined range of colors, to match up with the
# # number of species in the data 'tibl_podd' data frame
# clr04 <- clr03[1:length(unique(spcabbr.q))]
# col.f_spc <- clr04
# spcNm <- spcabbr.q
# # combine in to data frame
# df_clr <- as.data.frame(cbind(spcNm,col.f_spc))
# 
# # make vector with 3 species
# # as the first try out of the code is probably best to try out on
# # a limited number of species
# # spAbttry <- c("Hemsan",
# #               "Mnelei",
# #               "Psefar") 
# 
# 
# 
# 
# # Load libraries
# 
# library(ggplot2)
# library(dplyr)
# #install.packages("rstan")
# library(rstan)
# #install.packages("here")
# library(here);options(mc.cores = parallel::detectCores())
# 
# # Load data
# 
# getwd()
# 
# # make an object with the path to the data
# pth_gh <-        "https://raw.githubusercontent.com/gledguri/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples/d8c6b248cf8b0c57127b9005cc47ec686099fc6c/Data/"
# pf01 <- paste0(pth_gh,"/ddPCR_environmental_samples.csv")
# # list all files on the github data webpage
# f01 <- "Sample_metadata.csv"
# f02 <- "ddPCR_environmental_samples.csv"
# f03 <- "ddPCR_standard_samples.csv"
# f04 <- "ddpcr_amplitude_cod_herring_standard_samples.csv"
# f05 <- "ddpcr_amplitude_cod_saithe_standard_samples.csv"
# f06 <- "qPCR_environmental_samples.csv"
# f07 <- "qPCR_opt_Cod_Herring.csv"
# f08 <- "qPCR_opt_Cod_Saithe.csv"
# f09 <- "qPCR_standard_samples.csv"
# f10 <- "sample_idx.csv"
# # make a list of the files
# lrf <-list(f01,f02,f03,f04,f05,f06,f07,f08,f09,f10)
# # iterate over the list of files and read them into the environment
# for (i in lrf){
#   assign(i,read.csv(paste0(pth_gh,i)))
# }
# 
# 
# 
# 
# 
# # install.packages("curl")
# library(curl)
# # install.packages("RCurl")
# library(RCurl)
# 
# # 
# fls01 <- "R_code.qmd"
# fls02 <- "Stan_qPCR_ddPCR_prob_of_det.stan"
# fls03 <- "Stan_qPCR_ddPCR_quantification_precision.stan"
# fls04 <- "Stan_qPCR_two_step_model.stan"
# # make a list of the stan files from the github repository
# lst_stfls <- list(fls01,fls02,fls03,fls04)
# # iterate over the list of files and read them into the environment
# for (i in lst_stfls){
#   url_k <- paste0(pth_to_gh_stc,"/",i)
#   dest_pth_file <- paste0(outdir01,"/",i)
#   download.file(url_k,dest_pth_file,method="auto")
#   #readLines() %>% 
#   #writeLines(con = paste0(outdir01,"/",i))
# }
# 
# # Load data
# #install.packages("Rtools")
# #install.packages("rstan")
# library(rstan)
# library(here)
# 
# # Or try subsetting the df_ddpcr03 data frame to only comprise the
# # first 6 species 
# sbs_spcs <- spcabbr.q[1:2]
# # sbs_spcs <- spcabbr.q[1:5]
# # sbs_spcs <- spcabbr.q[1:6]
# 
# no_of_spsc <- seq(1,length(spcabbr.q),1)
# # split the number of species in sets of 3
# # to be able to iterate over 3 species per iteration
# noos <- split(no_of_spsc, ceiling(seq_along(no_of_spsc)/3))
# # get the number of elements in the list
# nsets <-length(noos)
# # make a sequence for the number of sets
# sqnsets <- seq(1,nsets,1)
# 
# #df_ddpcr03
# #df_qpcr03
# 
# # iterate over the list of set of species
# for (ns in sqnsets) {
#   print(ns)
# }
# 
# ns <- 4
# ns <- 5
# ns <- 6
#    #}
# # get the index numbers for the species
# nss <- unlist(noos[ns])
# #  use the index number for the species in the list of species
# # to get the species name
# sbs_spcs <- spcabbr.q[nss]
# sbs_spcs <- sbs_spcs[1]
# sbs_spcs
# #pad with zeros to two characters
# #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
# no_sbs <- stringr::str_pad(ns, 2, pad = "0")
# # match the  species names to get the colors that
# # correspond to the species
# clr04 <- df_clr$col.f_spc[match(sbs_spcs,df_clr$spcNm)]
# # get the species names for the subsetted species,
# # and concatenate into a string
# sbs_spcs_cnc <- paste0(sbs_spcs,collapse = "_")
# sbs_spcs_cnc <- sbs_spcs
# 
# # There are 17 species
# # I suspect that the Stanmod part will not work with 17 species
# 
# # collapse the vector to a string
# #testspc.Abbr <- paste0(spAbttry,collapse = "|")
# testspc.Abbr <- paste0(sbs_spcs,collapse = "|")
# 
# # subset the data frame to only comprise 
# # the species that are in the 'spAbttry' vector
# df_ddpcr04 <- df_ddpcr03[grepl(testspc.Abbr,df_ddpcr03$spcabbr),]
# # also subset the qpcr data frame to only comprise the
# # species that are in the 'spAbttry' vector
# df_qpcr04 <- df_qpcr03[grepl(testspc.Abbr,df_qpcr03$spcabbr),]
# 
# length(unique(df_qpcr04$spcabbr))
# length(unique(df_ddpcr04$spcabbr))
# # subset the ddpcr data frame to only comprise the 
# # the std and ntc samples, to have a data frame that 
# # only has the standard 
# # dilution series
# df_ddpcr.std <- df_ddpcr04[grepl("Standard",df_ddpcr04$WellType),]
# # and subset to only comprise the samples that are 
# # eDNA samples 
# df_ddpcr.e <- df_ddpcr04[grepl("Unknown",df_ddpcr04$WellType),]
# 
# # subset the qpcr data frame to only comprise the 
# # the std and ntc samples, to have a data frame 
# # that only has the standard 
# # dilution series
# df_qpcr.std <- df_qpcr04[grepl("Standard",df_qpcr04$WellType),]
# # and subset to only comprise the samples that are 
# # eDNA samples 
# df_qpcr.e <- df_qpcr04[grepl("Unknown",df_qpcr04$WellType),]
# 
# # check the number of unique species in the 'spcabbr' column
# # in the ddpcr data frame
# spa.in.dd <- unique(df_ddpcr.std$spcabbr)
# length(unique(spa.in.dd))
# spa.in.dd <- spa.in.dd[order(spa.in.dd)]
# # also check in the qpcr data frame
# spa.in.q <- unique(df_qpcr.std$spcabbr)
# length(unique(spa.in.q))
# spa.in.q <- spa.in.q[order(spa.in.q)]
# # check if the species are the same in the two data frames
# spa.in.q %in% spa.in.dd
# # identify the species that is not in both data frames
# spa.in.q[!spa.in.q %in% spa.in.dd]
# # it appears 'Hemtak' is missing from the ddpcr data frame !?
# # or that is , it was missing - The change of 'STD' to 'std'
# # made the difference, and made it possible to 
# # identify the 'Standard' samples for the WellType column
# 
# 
# # check what the different values for the std levels 
# # in the MST_sample_number column
# # looks like
# df_qpcr.std %>% dplyr::distinct(MST_sample_number) %>%
#   dplyr::arrange(MST_sample_number)
# df_ddpcr.std %>% dplyr::distinct(MST_sample_number) %>%
#   dplyr::arrange(MST_sample_number)
# 
# # The qPCR std data frame needs columns with
# # st_idx, species_idx and sp_st_idx, add this to the
# # 'std' data frame
# # add std_idx to the standard samples on the df_qpcr.std data frame
# df_qpcr.std <- df_qpcr.std %>%                                        # Create ID by group
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(st_idx = cur_group_id())
# # add species_idx to the standard samples on the df_qpcr.std data frame
# df_qpcr.std <- df_qpcr.std %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr) %>%
#   dplyr::mutate(species_idx = cur_group_id())
# # check the number of unique species in the 'species_idx' column
# length(unique(df_qpcr.std$species_idx))
# 
# # add a 'sp_st_idx' number to the standard samples on the q_ddpcr.std data frame
# df_qpcr.std <- df_qpcr.std %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr,st_idx) %>%
#   dplyr::mutate(sp_st_idx = cur_group_id())
# # check the unique st_idx, species_idx and sp_st_idx
# unique(df_qpcr.std$st_idx)
# unique(df_qpcr.std$species_idx)
# unique(df_qpcr.std$sp_st_idx)
# 
# # add a presence column to the standard samples on the
# # df_qpcr.std data frame
# # the presence column is a binary column that  is evaluated from the
# # 'Conc' column. The Conc column holds values of concentrations
# # translated from the Ct values. If the 'Conc' is above 0
# # then there is a presence, and the 'pres' column is evaluted to be '1'
# df_qpcr.std <- df_qpcr.std %>%                                        
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(pres = ifelse(Conc>0,1,0))
# 
# # The ddPCR std data frame also needs columns with
# # st_idx, species_idx and sp_st_idx, add this to the
# # 'std' data frame
# # add std_idx to the standard samples on the
# # df_ddpcr.std data frame
# df_ddpcr.std <- df_ddpcr.std %>%                                        # Create ID by group
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(st_idx = cur_group_id())
# # add species_idx to the standard samples on the df_ddpcr.std data frame
# df_ddpcr.std <- df_ddpcr.std %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr) %>%
#   dplyr::mutate(species_idx = cur_group_id())
# # add a 'sp_st_idx' number to the standard samples on the df_ddpcr.std data frame
# df_ddpcr.std <- df_ddpcr.std %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr,st_idx) %>%
#   dplyr::mutate(sp_st_idx = cur_group_id())
# # check the number of unique species in the 'species_idx' column
# length(unique(df_ddpcr.std$species_idx))
# 
# # add a presence column to the standard samples on
# # the df_qpcr.std data frame
# df_ddpcr.std <- df_ddpcr.std %>%                                        
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(pres = ifelse(Conc>0,1,0))
# 
# # The qPCR 'e' data frame also needs columns with 
# # 'st_idx', 'species_idx' and 'sp_st_idx', and 'pres'
# # add these columns to the df_qpcr.e data frame
# # using the approach as above
# # add std_idx to the standard samples on the df_qpcr.e data frame
# df_qpcr.e <- df_qpcr.e %>%                                        # Create ID by group
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(st_idx = cur_group_id())
# # add species_idx to the standard samples on the df_qpcr.e data frame
# df_qpcr.e <- df_qpcr.e %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr) %>%
#   dplyr::mutate(species_idx = cur_group_id())
# # check the number of unique species in the 'species_idx' column
# length(unique(df_qpcr.e$species_idx))
# # add a 'sp_st_idx' number to the standard samples on the df_qpcr.e data frame
# df_qpcr.e <- df_qpcr.e %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr,st_idx) %>%
#   dplyr::mutate(sp_st_idx = cur_group_id())
# # add a presence column to the standard samples 
# # on the df_qpcr.e data frame
# df_qpcr.e <- df_qpcr.e %>%                                        
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(pres = ifelse(Conc>0,1,0))
# 
# # The ddPCR 'e' data frame also needs columns with 
# # 'st_idx', 'species_idx' and 'sp_st_idx', and 'pres'
# # add these columns to the df_ddpcr.e data frame
# # using the approach as above
# # add std_idx to the standard samples on the df_ddpcr.e data frame
# df_ddpcr.e <- df_ddpcr.e %>%                                        # Create ID by group
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(st_idx = cur_group_id())
# # add species_idx to the standard samples on the df_ddpcr.e data frame
# df_ddpcr.e <- df_ddpcr.e %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr) %>%
#   dplyr::mutate(species_idx = cur_group_id())
# # check the number of unique species in the 'species_idx' column
# length(unique(df_ddpcr.e$species_idx))
# # add a 'sp_st_idx' number to the standard samples on the df_ddpcr.e data frame
# df_ddpcr.e <- df_ddpcr.e %>%                                        # Create ID by group
#   dplyr::group_by(spcabbr,st_idx) %>%
#   dplyr::mutate(sp_st_idx = cur_group_id())
# # add a presence column to the standard samples on the df_ddpcr.e data frame
# df_ddpcr.e <- df_ddpcr.e %>%                                        
#   dplyr::group_by(MST_sample_number) %>%
#   dplyr::mutate(pres = ifelse(Conc>0,1,0))
# 
# # check the unique values in the 'pres' column
# unique(df_ddpcr.e$pres)
# 
# 
# 
# 
# # the 4 data frames :
# # 'df_ddpcr.std'
# # 'df_ddpcr.e'
# # 'df_qpcr.std'
# # 'df_qpcr.e'
# 
# # should now equal the 
# # e_ddpcr <- ddPCR_environmental_samples.csv
# # e_qpcr <- qPCR_environmental_samples.csv
# # st_ddpcr <- ddPCR_standard_samples.csv
# # st_qpcr <- qPCR_standard_samples.csv
# # provided in the original code developed by Gledis Guri
# 
# 
# # use this string to subset the data frames
# df_ddpcr.std <- df_ddpcr.std
# df_ddpcr.e <- df_ddpcr.e
# df_qpcr.std <- df_qpcr.std
# df_qpcr.e <- df_qpcr.e
# 
# 
# st_ddpcr <- df_ddpcr.std
# st_qpcr <- df_qpcr.std
# # View(st_qpcr2)
# #e_ddpcr
# # View(st_ddpcr2)
# write.csv(df_ddpcr.std,"df_ddpcr.std.csv")
# write.csv(df_ddpcr.e,"df_ddpcr.e.csv")
# write.csv(df_qpcr.std,"df_qpcr.std.csv")
# write.csv(df_qpcr.e,"df_qpcr.e.csv")
# 
# #_____________________________________________________________________
# #_____________________________________________________________________
# #_____________________________________________________________________
# 
# # Quantifying the detection sensitivity and precision of qPCR and ddPCR 
# # mechanisms for eDNA samples
# # Author
# # 
# # Gled Guri
# # Published
# # from the scientific paper:
# # Guri, G., Ray.J.L., Shelton, A.O., Kelly, R.P., Præbel, K.,Allan, E.A., Yoccoz, N., Johansen, T., Wangensteen, O.S., Hanebrekke, T., Westgaard, J.-I., 2024. Quantifying the detection sensitivity and precision of qPCR and ddPCR mechanisms for eDNA samples. Ecology and Evolution, 14:e70678. https://doi.org/10.1002/ece3.70678
# # 
# # November 10, 2024
# # https://html-preview.github.io/?url=https://github.com/gledguri/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples/blob/main/Code/Quantifying-the-detection-sensitivity-and-precision-of-qPCR-and-ddPCR-mechanisms-for-eDNA-samples.html
# 
# # Load functions
# 
# extract_list_param <- function(stanMod){
#   l <- stanMod@model_pars
#   x <- list(extract_param(stanMod,"lp__")) %>% setNames("lp__")
#   for (i in l) {
#     x <- c(x,list(extract_param(stanMod,i)) %>% setNames(i))
#   }
#   x <- x[-which(names(x) == c("lp__"))]
#   return(x)
# }
# 
# cloglog <- function(theta) log(-log(1 - theta))
# 
# extract_param <- function(model=stanMod,parmeter="alpha"){
#   return(summary(model, par = parmeter)$summary %>% 
#            unlist()%>%as.data.frame%>%round(.,2))
# }
# 
# inv.cloglog <- function(theta) 1-exp(-exp(theta))
# 
# logreg <- function(x,b0,b1,i) y=(1/(1+exp(-(b0+(b1*(x+i))))))
# 
# scientific_10 <- function(x) {
#   c <- scales::scientific_format()(x)
#   t <- gsub("1e", "10^", c)
#   t2 <- gsub("10\\^\\+", "10\\^", t)
#   str2expression(t2)}
# 
# # Load data
# 
# # e_ddpcr <- read.csv(here('Data','ddPCR_environmental_samples.csv'))
# # e_qpcr <- read.csv(here('Data','qPCR_environmental_samples.csv'))
# # st_ddpcr <- read.csv(here('Data','ddPCR_standard_samples.csv'))
# # st_qpcr <- read.csv(here('Data','qPCR_standard_samples.csv'))
# 
# # e_ddpcr <- ddPCR_environmental_samples.csv
# # e_qpcr <- qPCR_environmental_samples.csv
# # st_ddpcr <- ddPCR_standard_samples.csv
# # st_qpcr <- qPCR_standard_samples.csv
# 
# st_ddpcr <- df_ddpcr.std
# st_qpcr <- df_qpcr.std
# e_ddpcr <- df_ddpcr.e
# e_qpcr <- df_qpcr.e
# # use dplyr filter to remove the rows that have NAs in the 'Ct' column
# st_qpcr_cm <- st_qpcr %>% dplyr::filter(!is.na(Ct))
# e_qpcr_cm <- e_qpcr %>% dplyr::filter(!is.na(Ct))
# 
# # sample_metadata <- read.csv(here('Data','Sample_metadata.csv'))
# # stidx <- read.csv(here('Data','sample_idx.csv'))
# 
# # sample_metadata <- Sample_metadata.csv
# # stidx <- sample_idx.csv
# 
# #View(st_qpcr)
# #View(e_qpcr)
# # 1 Probability of detection
# # 1.1 Create stan data
# 
# stan_data_M1 <- list(
#   N_i = length(unique(st_qpcr$Species)),
#   
#   Nq = nrow(st_qpcr),
#   Ndd = nrow(st_ddpcr),
#   
#   i_q_idx = st_qpcr$species_idx,
#   i_d_idx = st_ddpcr$species_idx,
#   
#   Z_qPCR = st_qpcr$pres,
#   Z_ddPCR = st_ddpcr$pres,
#   
#   # C_qPCR = log10(st_qpcr$int_concentation),
#   # C_ddPCR = log10(st_ddpcr$int_concentation)
#   C_qPCR = log10(st_qpcr$init_concentration),
#   C_ddPCR = log10(st_ddpcr$init_concentration)
# );str(stan_data_M1)
# 
# # 1.2 Create initial values
# initial_values <- replicate(4,list(
#   phi_qPCR_0 = rep(1.0, 3),
#   phi_qPCR_1 = rep(0.4, 3),
#   phi_ddPCR_0 = rep(2.50, 3),
#   phi_ddPCR_1 = rep(1.23, 3)
# ),simplify=F)
# 
# initial_values
# # 1.2 The math
# 
# ## The stan code
# # 
# # This stan code runs 2 independent logistic regression analyses
# # (equation 1 and 2) of qPCR and ddPCR data, with separate sets of 
# # parameters 
# # (intercepts and slopes) for each technique. 
# # To ensure comparability between the parameters
# # , and , , we adjusted by adding 2 in stan ( copies/μL) since
# # the concentration
# # of ddPCR standards ranged from to , while those for qPCR 
# # range from to .
# # The shift aligns the slope and intercept enabling to 
# # compare the two models 
# # through their parameters ( and).
# 
# # #______________________________________________________________________
# # # I am unable to get the code below to work. 
# # #______________________________________________________________________
# # # Gled: This part is the Stan language (C++) code for running the Bayesian Model. 
# # # This can't be run in R but it can be called from R and run through rstan package (line 198-205). Skip lines 151 - 194 an move straight to 1.3 Run stan Model
# # # The code below is in 'Stan_qPCR_ddPCR_prob_of_det.stan' and it will be called it in line 199. 
# # # I left it in markdown format for educational purpose only.
# # data { 
# #   //Numbers of dimentions
# #   int Nq; // Total number of observation in qPCR standard samples
# #   int Ndd; // Total number of observation in ddPCR standard samples
# #   int N_i; // Number of assays in both data (qPCR and ddPCR)
# #   //Indexes   
# #   array[Nq] int i_q_idx; // Species/assay index for qPCR 
# #   array[Ndd] int i_d_idx; // Species/assay index for ddPCR
# #   // Data
# #   array[Nq] int Z_qPCR; // Presence/Absence of targets in qPCR runs
# #   array[Ndd] int Z_ddPCR; // Presence/Absence of targets in ddPCR runs
# #   array[Nq] real C_qPCR; // Known concentration (log10) in qPCR data
# #   array[Ndd] real C_ddPCR; // Known concentration (log10) in qPCR data
# # }
# # parameters {
# #   vector[N_i] phi_qPCR_0; // Intercept of logistic regression for qPCR detection probability model
# #   vector[N_i] phi_qPCR_1; // Slope of logistic regression for qPCR detection probability model
# #   vector[N_i] phi_ddPCR_0; // Intercept of logistic regression for ddPCR detection probability model
# #   vector[N_i] phi_ddPCR_1; // Slope of logistic regression for ddPCR detection probability model
# # }
# # transformed parameters{
# #   vector[Nq] theta_qPCR; // probability of detection (logit-space) qPCR model
# #   vector[Ndd] theta_ddPCR; // probability of detection (logit-space) ddPCR model
# #   // Logistic regression for qPCR detection probability model
# #   for (i in 1:Nq){
# #     theta_qPCR[i] = phi_qPCR_0[i_q_idx[i]] + (phi_qPCR_1[i_q_idx[i]] * C_qPCR[i]);
# #   }
# #   // Logistic regression for ddPCR detection probability model
# #   for (i in 1:Ndd){
# #     theta_ddPCR[i] = phi_ddPCR_0[i_d_idx[i]] + (phi_ddPCR_1[i_d_idx[i]] * (C_ddPCR[i]+2));
# #     #//' Since the standards of ddPCR range from 10^-3 - 10^4 while qPCR 
# #     #//' standards range from 10^-1 - 10^6 we added 2 to C_ddPCR to make the 
# #     #//' parameters phi_qPCR_0 and phi_qPCR_1 comparable with phi_ddPCR_0 and phi_ddPCR_1
# #   }
# # }
# # model {
# #   Z_qPCR ~ bernoulli(inv_logit(theta_qPCR)); 
# #   Z_ddPCR ~ bernoulli(inv_logit(theta_ddPCR));
# #   // Priors
# #   phi_qPCR_0 ~ normal(0, 1);
# #   phi_qPCR_1 ~ normal(0, 1);
# #   phi_ddPCR_0 ~ normal(0, 1);
# #   phi_ddPCR_1 ~ normal(0, 1);
# # }
# # #
# # 
# # #______________________________________________________________________
# # #______________________________________________________________________
# 
# # 1.3 Run stan Model
# 
# stanMod_M1 = stan(file = here::here('code',
#                                     'Stan_qPCR_ddPCR_prob_of_det.stan'),
#                   chains = 4,
#                   model_name = 'qPCR_vs_ddPCR_probability_of_detection',
#                   iter = 5000,
#                   warmup = 2000,
#                   # init = initial_values,
#                   data = stan_data_M1)
# 
# #1.4 Data for Figure 1
# 
# # Extract parameters from stan model to build the graph
# mod_out <- extract_list_param(stanMod_M1)
# mod_out <- mod_out[which(names(mod_out) %in%
#                            c("phi_qPCR_0",
#                              "phi_qPCR_1",
#                              "phi_ddPCR_0",
#                              "phi_ddPCR_1"))]
# # make a sequence of the species index
# sp_idx.sq <- (unique(st_ddpcr$species_idx))
# length(sp_idx.sq)
# sp_idx.no <- (unique(st_qpcr$species_idx))
# spcabr.nm <- (unique(st_qpcr$spcabbr))
# df_sp.sq <- as.data.frame(cbind(sp_idx.no,spcabr.nm))
# df_sp.sq$sp_idx.no <- as.numeric(df_sp.sq$sp_idx.no)
# # count the number of species
# # this count will be used later on for producing the
# # probability_of_detection_qpcr and the
# # probability_of_detection_ddpcr data frames
# Nspc <- length(sp_idx.no)
# # reorder the data frame by the  index number for the
# # species
# df_sp.sq <- df_sp.sq[order(df_sp.sq$sp_idx.no),]
# # us this data frame for the 'mod_out' loop
# # Add species names
# for (i in names(mod_out)) {
#   mod_out[[i]] <- mod_out[[i]] %>%
#     cbind(.,tibble(#sp_idx = 1:3,
#       sp_idx =c(df_sp.sq$sp_idx.no),
#                    #Species = c("Cod", "Herring", "Saithe")
#       Species = c(df_sp.sq$spcabr.nm)
#     ))
# }
# 
# # Construct the logistic regression for qPCR probability of 
# # detection from phi_qPCR_0 & phi_qPCR_1 parameters
# probability_of_detection_qpcr <-
#   as.data.frame(matrix(rep(seq(-3,6,by=0.1),Nspc),
#                        ncol=1)) %>% setNames("C") %>%
#   dplyr::mutate(Species=rep(
#     c(df_sp.sq$spcabr.nm),
#     #c("Cod","Herring","Saithe"),
#                             each=nrow(.)/Nspc)) %>%
#   left_join(.,mod_out$phi_qPCR_0 %>% 
#               select(Species,mean),by="Species") %>%
#   dplyr::rename(intercept=mean) %>%
#   dplyr::left_join(.,mod_out$phi_qPCR_1 %>% 
#                      select(Species,mean),by="Species") %>%
#   dplyr::rename(slope=mean) %>%
#   dplyr::mutate(prob=logreg(C,intercept,slope,0))
# 
# # Construct the logistic regression for ddPCR probability of 
# # detection from alpha_0 & alpha_1 parameters
# probability_of_detection_ddpcr <-
#   as.data.frame(matrix(rep(seq(-3,6,by=0.1),Nspc),ncol=1)) %>% 
#   setNames("C") %>%
#   dplyr::mutate(Species=rep(
#     c(df_sp.sq$spcabr.nm),
#     #c("Cod","Herring","Saithe"),
#                             each=nrow(.)/Nspc)) %>%
#   dplyr::left_join(.,mod_out$phi_ddPCR_0 %>% 
#                      select(Species,mean),by="Species") %>%
#   dplyr::rename(intercept=mean) %>%
#   dplyr::left_join(.,mod_out$phi_ddPCR_1 %>% 
#                      select(Species,mean),by="Species") %>%
#   dplyr::rename(slope=mean) %>%
#   dplyr::mutate(prob=logreg(C,
#                             #added 2 to logistic regression (see why in stan model file)
#                             intercept,slope,2)) 
# 
# # Combine the logistic regresion lines (prob of det) of both data (qPCR and ddPCR)
# probability_of_detection_data <- cbind(
#   probability_of_detection_qpcr %>% 
#     dplyr::select(C,Species,prob) %>%
#     dplyr::rename(prob_q=prob),
#   probability_of_detection_ddpcr %>% 
#     dplyr::select(prob) %>% 
#     dplyr::rename(prob_dd=prob)) %>%
#   dplyr::mutate(delta=prob_dd-prob_q)
# 
# probability_of_detection_data %>% 
#   as_tibble() %>% 
#   print(n=20)
# tibl_podd <- probability_of_detection_data
# 
# # # A tibble: 273 × 5
# #        C Species  prob_q prob_dd delta
# #    <dbl> <chr>     <dbl>   <dbl> <dbl>
# #  1  -3   Cod     0.00875   0.209 0.200
# #  2  -2.9 Cod     0.0106    0.233 0.222
# #  3  -2.8 Cod     0.0128    0.258 0.246
# #  4  -2.7 Cod     0.0155    0.286 0.270
# #  5  -2.6 Cod     0.0187    0.315 0.296
# #  6  -2.5 Cod     0.0225    0.345 0.323
# #  7  -2.4 Cod     0.0272    0.377 0.350
# #  8  -2.3 Cod     0.0327    0.410 0.377
# #  9  -2.2 Cod     0.0394    0.444 0.404
# # 10  -2.1 Cod     0.0473    0.478 0.431
# # 11  -2   Cod     0.0568    0.512 0.456
# # 12  -1.9 Cod     0.0680    0.547 0.479
# # 13  -1.8 Cod     0.0812    0.581 0.500
# # 14  -1.7 Cod     0.0967    0.614 0.517
# # 15  -1.6 Cod     0.115     0.646 0.531
# # 16  -1.5 Cod     0.136     0.677 0.541
# # 17  -1.4 Cod     0.160     0.706 0.546
# # 18  -1.3 Cod     0.188     0.734 0.547
# # 19  -1.2 Cod     0.219     0.760 0.542
# # 20  -1.1 Cod     0.253     0.784 0.531
# # ℹ 253 more rows
# 
# 
# # 1.5 Plot Figure 1 - Detection probability (qPCR vs ddPCR)
# 
# # source(here('code','Surpressed','Figure_1.R'))
# # Figure_1
# # 
# # # get one range of colors
# # clr01 <- palette.colors(palette = "Okabe-Ito")
# # # get another range of colors
# # clr02 <- palette.colors(palette = "Polychrome")
# # # combine the two ranges of colors
# # clr03 <- c(clr01,clr02)
# # # use the combined range of colors, to match up with the
# # # number of species in the data 'tibl_podd' data frame
# # clr04 <- clr03[1:length(unique(tibl_podd$Species))]
# 
# #pad with zeros to two characters
# #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
# nos <- stringr::str_pad(ns, 2, pad = "0")
# # concatenate the species names, to be able to use them
# # for creating a file name
# sbs_spcs_cnc <- paste0(sbs_spcs,collapse = "_")
# # save the tibble with the detection probability data
# # in the wd00_14 directory as a csv file
# # with the concatenated species names
# # and the number for the subset of species
# write.csv(tibl_podd,paste0(outdir02,
#                            "/tibl_podd_",no_sbs,"_",
#                            sbs_spcs_cnc,".csv"))
# 
# 
# # make a ggplot object
# pp1 <-
#   tibl_podd %>%
#   ggplot()+
#   geom_line(aes(x=10^C,y=prob_dd,color=Species),lty=2,size=0.9)+
#   geom_line(aes(x=10^C,y=prob_q,color=Species),lty=1,size=0.9)+
#   # geom_point(data=st_qpcr %>% group_by(Species,Sample_name) %>%
#   #              summarise(int_concentation=mean(int_concentation),
#   #                        pres=mean(pres)),
#   #            aes(x=int_concentation,y=pres,color=Species),shape=15,size=3)+
#   # geom_point(data=st_ddpcr %>% group_by(Species,Sample_name) %>%
#   #              summarise(int_concentation=mean(int_concentation),
#   #                        pres=mean(pres)),
#   #            aes(x=int_concentation,y=pres,color=Species),shape=19,size=3)+
#   # scale_x_log10(labels=scientific_10,lim=c(1e-4,1e6))+
#   scale_x_log10(labels=NULL,lim=c(1e-2,1e6),
#                 breaks=10^c(-2,-1,0,1,2,4,6))+
#   ylim(0,1)+
#   scale_color_manual(values=
#                                     
#                        #c("tomato2", "deepskyblue2","orange2"))+
#   # in the original code only 3 species were included
#   # this only required 3 colors. This setup has 17 species
#   # and need 17 colors
#                        c(clr04))+
#     
#   facet_grid(~Species)+
#   theme_bw()+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         #strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         #plot.margin = margin(0.1, 0.1, -0.5, 0.5, "cm"),
#         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
#         axis.title.x = element_text(size = 19),
#         axis.text.x = element_text(size=15),
#         axis.text.y=element_text(size=15))+
#   ylab("Detection probability")+
#   xlab("C - Nominal concentration (copies/μL)")
# #xlab("")
# pp1
# 
# #set variable to define if figures are to be saved
# bSaveFigures<-T
# # save the figure if the above 'bSaveFigures' is TRUE
# if(bSaveFigures==T){
#   ggsave(plot = pp1, 
#          # define the output filenmae by pasting together 
#          # qpcrrunno and qpcrrundate
#          filename = paste0(outdir02,
#                            "/Fig01_1_Detection_probability_",
#                            no_sbs,"_",sbs_spcs_cnc,
#                            ".png"),
#          width=210*1.6,height=297*0.6,
#          #width=297,height=210,
#          units="mm",dpi=300)
# }
# 
# pp1 <-
#   tibl_podd %>%
#   ggplot()+
#   geom_line(aes(x=10^C,y=prob_dd,color=Species),lty=2,size=0.9)+
#   geom_line(aes(x=10^C,y=prob_q,color=Species),lty=1,size=0.9)+
#   # geom_point(data=st_qpcr %>% group_by(Species,Sample_name) %>%
#   #              summarise(int_concentation=mean(int_concentation),
#   #                        pres=mean(pres)),
#   #            aes(x=int_concentation,y=pres,color=Species),shape=15,size=3)+
#   # geom_point(data=st_ddpcr %>% group_by(Species,Sample_name) %>%
#   #              summarise(int_concentation=mean(int_concentation),
#   #                        pres=mean(pres)),
#   #            aes(x=int_concentation,y=pres,color=Species),shape=19,size=3)+
#   # scale_x_log10(labels=scientific_10,lim=c(1e-4,1e6))+
#   scale_x_log10(labels=NULL,lim=c(1e-2,1e6),breaks=10^c(-2,-1,0,1,2,4,6))+
#   ylim(0,1)+
#   #scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
#   scale_color_manual(values=c(clr04))+
#   facet_grid(~Species)+
#   theme_bw()+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         plot.margin = margin(0.1, 0.1, -0.5, 0.5, "cm"),
#         #plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
#         axis.title.x = element_text(size = 19),
#         axis.text.x = element_text(size=15),
#         axis.text.y=element_text(size=15))+
#   ylab("Detection probability")+
#   #xlab("C - Nominal concentration (copies/μL)")
#   xlab("")
# pp1
# 
# theta_line <- tibl_podd
# pp2 <-
#   theta_line %>%
#   ggplot()+
#   geom_line(aes(x=10^C,y=delta,color=Species),
#             lty=6,size=0.9)+
#   scale_x_log10(labels=scientific_10,
#                 lim=c(1e-2,1e6),
#                 breaks=10^c(-2,-1,0,1,2,4,6))+
#   ylim(-0.1,0.8)+
#   #scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
#   scale_color_manual(values=c(clr04))+
#   facet_grid(~Species)+
#   theme_bw()+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
#         axis.title.x = element_text(size = 19),
#         axis.text.x = element_text(size=15),
#         axis.text.y=element_text(size=15))+
#   ylab("ddPCR - qPCR \ndetection probability")+
#   xlab("Nominal DNA concentration (copies/uL)")
# 
# pp2
# 
# # save the figure if the above 'bSaveFigures' is TRUE
# if(bSaveFigures==T){
#   ggsave(plot = pp2, 
#          # define the output filenmae by pasting together 
#          filename = paste0(outdir02,
#                            "/Fig01_2_Nominal_DNA_concentration_",
#                            no_sbs,"_",sbs_spcs_cnc,
#                            ".png"),
#          width=210*1.6,height=297*0.6,
#          #width=297,height=210,
#          units="mm",dpi=300)
# }
# 
# p1_leg <- ggplot() +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[1]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[2]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[3]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[4]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[5]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[6]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[7]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[8]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[9]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[10]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[11]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[12]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[13]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[14]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[15]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[16]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[17]), size = 1) +
#   
#   #
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[1]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[2]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[3]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[4]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[5]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[6]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[7]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[8]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[9]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[10]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[11]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[12]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[13]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[14]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[15]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[16]), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[17]), size = 1) +
#   #
#   theme_bw()+
#   scale_color_manual(
#     name='Assay',
#     #breaks=c("Gadus morhua", "Clupea harengus","Pollachius virens"),
#     breaks=c(df_sp.sq$spcabr.nm),
#     #values = c("tomato2", "deepskyblue2","orange"),
#     values = c(clr04),
#     guide = guide_legend(
#       # override.aes = list(lty = c(1,1,1),
#       #                     size = c(3,3,3)))
#     override.aes = list(lty = c(rep(1,Nspc) ),
#                         size = c(rep(3,Nspc))))
#     ) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16,face = "italic"))
# 
# 
# p1_leg_leg <- cowplot::get_legend(p1_leg+
#                         theme(legend.justification = 
#                         c(0,-1.3)))
# 
# 
# p1_leg_leg
# p3_leg <- ggplot() +
#   geom_line(aes(x = NA, y = NA, color = "ddPCR sensitivity"), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = "qPCR sensitivity"), size = 1) +
#   theme_bw()+
#   scale_color_manual(
#     name="Modelled sensitivity",
#     breaks=c("ddPCR sensitivity","qPCR sensitivity"),
#     values = c("black","black"),
#     guide = guide_legend(
#       override.aes = list(lty = c(2,1),
#                           size = c(1,2)))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16))
# p3_leg_leg <- cowplot::get_legend(p3_leg+
#                                     theme(legend.justification = 
#                                             c(0,-0.1)))
# 
# p2_leg <- ggplot() +
#   geom_line(aes(x = NA, y = NA, 
#                 color = "Sensitivity difference\n(ddPCR - qPCR)"), 
#             size = 1) +
#   theme_bw()+
#   scale_color_manual(
#     name="",
#     breaks=c("Sensitivity difference\n(ddPCR - qPCR)"),
#     values = c("black"),
#     guide = guide_legend(
#       override.aes = list(lty = c(6),
#                           size = c(1)))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=1),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16))
# p2_leg_leg <- cowplot::get_legend(p2_leg+
#       # the 'theme(legend.justification' appears
#         # to be able to shift the legends on side more around
#         # it is here changed from 0.95 to 10.95
#                                     theme(legend.justification = 
#                                             c(0,10.95)))
# 
# p4_leg <-
#   ggplot() +
#   geom_point(aes(x = NA, y = NA, color = "ddPCR standard samples"), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = "qPCR standard samples"), size = 1) +
#   theme_bw()+
#   scale_color_manual(
#     name="Measured sensitivity",
#     breaks=c("ddPCR standard samples","qPCR standard samples"),
#     values = c("black","black"),
#     guide = guide_legend(
#       override.aes = list(shape = c(19,15),
#                           size = c(3)))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16))
# p4_leg_leg <- cowplot::get_legend(p4_leg+
#                                     theme(legend.justification = 
#                                             c(0,2.2)))
# 
# legend <- cowplot::plot_grid(p1_leg_leg,
#                              p3_leg_leg,
#                              p2_leg_leg,
#                              p4_leg_leg,
#                              nrow = 4)
# 
# p1 <- cowplot::plot_grid(pp1,pp2,nrow = 2,
#                          rel_heights = c(4.7,5.3),
#                          labels = c("a","b"),
#                          label_size = 20, 
#                          align = "v")
# Figure_1 <- cowplot::plot_grid(p1,
#                                legend,nrow = 1,
#                                ncol = 2,
#                                rel_widths = c(7,1.8))
# Figure_1
# 
# # save the figure if the above 'bSaveFigures' is TRUE
# if(bSaveFigures==T){
#   ggsave(plot = Figure_1, 
#          # define the output filenmae by pasting together 
#          filename = paste0(outdir02,
#             "/Fig01_03_Nominal_DNA_",
#             "concentration_and_detected_for_",
#                            no_sbs,"_",sbs_spcs_cnc,
#                            ".png"),
#          width=210*2.4,height=297*1.2,
#          #width=297,height=210,
#          units="mm",dpi=300)
# }
# #_________________________________________________________________________________
# # 2 Quantification precision
# # 2.1 Create stan data
# 
# stan_data_M2 <- list(
#   N_st_q = nrow(st_qpcr),
#   N_en_q = nrow(e_qpcr),
#   N_st_qp = nrow(st_qpcr_cm),
#   N_en_qp = nrow(e_qpcr_cm),
#   #
#   N_st_d = nrow(st_ddpcr),
#   N_en_d = nrow(e_ddpcr),
#   #
#   N_i = length(unique(st_qpcr$Species)),
#   N_ij = length(unique(e_ddpcr$sp_st_idx)),
#   #
#   i_qst_idx = st_qpcr$species_idx,
#   i_qen_idx = e_qpcr$species_idx,
#   ij_qen_idx = e_qpcr$sp_st_idx,
#   
#   i_qst_p_idx = st_qpcr_cm$species_idx,
#   i_qen_p_idx = e_qpcr_cm$species_idx,
#   ij_qen_p_idx = e_qpcr_cm$sp_st_idx,
#   
#   Z_qPCR_st = st_qpcr$pres,
#   Z_qPCR_en = e_qpcr$pres,
#   #C_qPCR_st = log10(st_qpcr$int_concentation),
#   C_qPCR_st = log10(st_qpcr$init_concentration),
#   
#   Y_qPCR_st = st_qpcr_cm$Ct,
#   Y_qPCR_en = e_qpcr_cm$Ct,
#   #C_qPCR_st_continuous = log10(st_qpcr_cm$int_concentation),
#   C_qPCR_st_continuous = log10(st_qpcr_cm$init_concentration),
#   #
#   W_st = st_ddpcr$Positives,
#   W_en = e_ddpcr$Positives,
#   #
#   U_st = st_ddpcr$Tot_drop,
#   U_en = e_ddpcr$Tot_drop,
#   #
#   #C_ddPCR_st = log10(st_ddpcr$int_concentation),
#   C_ddPCR_st = log10(st_ddpcr$init_concentration),
#   #
#   i_dst_idx = st_ddpcr$species_idx,
#   #
#   i_den_idx = e_ddpcr$species_idx,
#   ij_den_idx = e_ddpcr$sp_st_idx
# );str(stan_data_M2)
# 
# # List of 28
# #  $ N_st_q              : int 277
# #  $ N_en_q              : int 1379
# #  $ N_st_qp             : int 242
# #  $ N_en_qp             : int 663
# #  $ N_st_d              : int 247
# #  $ N_en_d              : int 1052
# #  $ N_i                 : int 3
# #  $ N_ij                : int 343
# #  $ i_qst_idx           : int [1:277] 1 1 1 1 1 1 1 1 1 1 ...
# #  $ i_qen_idx           : int [1:1379] 1 1 1 1 2 2 2 2 1 1 ...
# #  $ ij_qen_idx          : int [1:1379] 1 1 1 1 172 172 172 172 2 2 ...
# #  $ i_qst_p_idx         : int [1:242] 1 1 1 1 1 1 1 1 1 1 ...
# #  $ i_qen_p_idx         : int [1:663] 1 1 1 2 2 1 2 2 2 2 ...
# #  $ ij_qen_p_idx        : int [1:663] 1 1 2 173 173 3 174 175 175 175 ...
# #  $ Z_qPCR_st           : int [1:277] 1 1 1 1 1 1 1 1 1 1 ...
# #  $ Z_qPCR_en           : int [1:1379] 1 1 0 0 0 0 0 0 0 0 ...
# #  $ C_qPCR_st           : num [1:277] 6 5 4 3 2 1 0 6 4 3 ...
# #  $ Y_qPCR_st           : num [1:242] 17.2 20.3 24.3 26.9 30.7 ...
# #  $ Y_qPCR_en           : num [1:663] 43.9 46.2 49.6 40.6 40.4 ...
# #  $ C_qPCR_st_continuous: num [1:242] 6 5 4 3 2 1 0 6 4 3 ...
# #  $ W_st                : int [1:247] 0 1 20070 20292 19709 19933 19559 19757 19520 19679 ...
# #  $ W_en                : int [1:1052] 0 0 0 0 0 0 0 0 3 2 ...
# #  $ U_st                : int [1:247] 20155 20155 20299 20299 19939 19939 19760 19760 19691 19691 ...
# #  $ U_en                : int [1:1052] 18263 18727 18377 18263 18727 18377 18711 18411 17584 18711 ...
# #  $ C_ddPCR_st          : num [1:247] -3 -3 4 4 4 4 4 4 4 4 ...
# #  $ i_dst_idx           : int [1:247] 2 1 2 1 2 1 2 1 2 1 ...
# #  $ i_den_idx           : int [1:1052] 1 1 1 2 2 2 1 1 1 2 ...
# #  $ ij_den_idx          : int [1:1052] 1 1 1 172 172 172 2 2 2 173 ...
# 
# initial_values <- replicate(4,list(
#   kappa_0 = rep(-11, Nspc),
#   kappa_1 = rep(0.48, Nspc)
# ),simplify=F)
# 
# #  2.2 The math
# 
# ## The stan code This stan code runs 2 independent models (equation 1 - 5 and 6 - 7) for qPCR and ddPCR data. Both models are run simultaneously for standard and environmental samples where the known concentration from standards informs the intercept and the slope parameters which thereafter informs the initial concentration of unknown samples (C).
# 
# # data {
# #   //Numbers of dimensions
# #   // // // qPCR
# #   int N_st_q; // Total number of observation in qPCR standard samples
# #   int N_en_q; // Total number of observation in qPCR environmental samples
# #   int N_st_qp; // Total number of observation in qPCR standard samples for only detected samples
# #   int N_en_qp; // Total number of observation in qPCR environmental samples for only detected samples
# #   // // // ddPCR
# #   int N_st_d; // Total number of observation in ddPCR standard samples
# #   int N_en_d; // Total number of observation in qPCR environmental samples
# #   // // Joined
# #   int N_i; // Number of species in both data
# #   int N_ij; // Number of species and stations in both data
# #   //
# #     //Indexes
# #   // // qPCR
# #   // // // Binomial model
# #   array[N_st_q] int i_qst_idx; // Species index for qPCR standard samples
# #   array[N_en_q] int i_qen_idx; // Species index for qPCR environmental samples
# #   array[N_en_q] int ij_qen_idx; // Species and standard index for qPCR environmental samples
# #   // // // Continious model
# #   array[N_st_qp] int i_qst_p_idx; // Species index for qPCR standard samples
# #   array[N_en_qp] int i_qen_p_idx; // Species index for qPCR environmental samples
# #   array[N_en_qp] int ij_qen_p_idx; // Species and standard index for qPCR environmental samples
# #   // // ddPCR
# #   array[N_st_d] int i_dst_idx; // Species index for ddPCR environmental samples
# #   array[N_en_d] int i_den_idx; // Species index for ddPCR standard samples
# #   array[N_en_d] int ij_den_idx; // Species and standard index for ddPCR environmental samples
# #   //
# #     // Data
# #   // // qPCR
# #   // // // Binomial model
# #   array[N_st_q] int Z_qPCR_st; // Presence/Absence response of qPCR standard data
# #   array[N_en_q] int Z_qPCR_en; // Presence/Absence response of qPCR environmental data
# #   array[N_st_q] real C_qPCR_st; // Known concentration (log10) in qPCR data
# #   // // // Continuous model
# #   array[N_st_qp] real Y_qPCR_st; // Ct values of qPCR standard data for only detected samples
# #   array[N_en_qp] real Y_qPCR_en; // Ct values of qPCR environmental data for only detected samples
# #   array[N_st_qp] real C_qPCR_st_continuous; // Known concentration (log10) in qPCR data for only detected samples
# #   // // ddPCR
# #   array[N_st_d] int W_st; // Observed positive droplets in ddPCR standard samples
# #   array[N_en_d] int W_en; // Observed positive droplets in ddPCR environmental samples
# #   array[N_st_d] int U_st; // Total droplets in ddPCR standard samples
# #   array[N_en_d] int U_en; // Total droplets in ddPCR environmental samples
# #   array[N_st_d] real C_ddPCR_st; // Known concentration (log10) in ddPCR data
# # }
# # parameters {
# #   // Parameters
# #   // // qPCR
# #   // // // Bernoulli model
# #   vector[N_i] phi_0;
# #   vector[N_i] phi_1;
# #   // // // Continous model
# #   vector[N_i] beta_0;
# #   vector[N_i] beta_1;
# #   vector[N_i] gamma_0;
# #   vector<upper=0>[N_i] gamma_1;
# #   vector[N_ij] C_qPCR;
# #   // // ddPCR
# #   vector[N_i] kappa_0;
# #   vector[N_i] kappa_1;
# #   vector<lower=-7>[N_ij] C_ddPCR;
# # }
# # transformed parameters{
# #   // Parameters
# #   // // qPCR
# #   // // // Bernoulli model
# #   vector[N_st_q] theta_st;
# #   vector[N_en_q] theta_un;
# #   // // // Continuous model
# #   vector[N_st_qp] mu_st;
# #   vector[N_en_qp] mu_en;
# #   vector[N_st_qp] sigma_st;
# #   vector[N_en_qp] sigma_en;
# #   // // ddPCR
# #   vector[N_st_d] omega_st;
# #   vector[N_en_d] omega_en;
# #   vector[N_ij] delta; // Difference between ddPCR and qPCR concentration estimates
# #   //
# #     // Model TP
# #   // // qPCR model
# #   // // // Bernuli module model compartment
# #   // // // // // Standard
# #   for (i in 1:N_st_q){
# #     theta_st[i] = phi_0[i_qst_idx[i]] + (phi_1[i_qst_idx[i]] * C_qPCR_st[i]);
# #   }
# #   // // // // // Unknown (Env samples)
# #   for (i in 1:N_en_q){
# #     theta_un[i] = phi_0[i_qen_idx[i]] + (phi_1[i_qen_idx[i]] * C_qPCR[ij_qen_idx[i]]);
# #   }
# #   // // // Continuous model compartment
# #   // // // // Standard
# #   for (i in 1:N_st_qp){
# #     mu_st[i] = beta_0[i_qst_p_idx[i]] + (beta_1[i_qst_p_idx[i]] * C_qPCR_st_continuous[i]);
# #     sigma_st[i] = exp(gamma_0[i_qst_p_idx[i]]+(gamma_1[i_qst_p_idx[i]] * C_qPCR_st_continuous[i]));
# #   }
# #   // // // // Unknown (Env samples)
# #   for (i in 1:N_en_qp){
# #     mu_en[i] = beta_0[i_qen_p_idx[i]] + (beta_1[i_qen_p_idx[i]] * C_qPCR[ij_qen_p_idx[i]]);
# #     sigma_en[i] = exp(gamma_0[i_qen_p_idx[i]]+(gamma_1[i_qen_p_idx[i]] * C_qPCR[ij_qen_p_idx[i]]));
# #   }
# #   // ddPCR Model
# #   // // Standard
# #   for (i in 1:N_st_d){
# #     omega_st[i] = kappa_0[i_dst_idx[i]]+(kappa_1[i_dst_idx[i]]*C_ddPCR_st[i]);
# #   }
# #   // // // Unknown (Env samples)
# #   for (i in 1:N_en_d){
# #     omega_en[i] = kappa_0[i_den_idx[i]]+(kappa_1[i_den_idx[i]]*C_ddPCR[ij_den_idx[i]]);
# #   }
# #   for (i in 1:N_ij){
# #     delta[i] = C_ddPCR[i]-C_qPCR[i];
# #   }
# # }
# # model {
# #   // Model
# #   // // qPCR
# #   // // // Bernoulli model
# #   Z_qPCR_st ~ bernoulli(inv_logit(theta_st)); //Standards
# #   Z_qPCR_en ~ bernoulli(inv_logit (theta_un)); //Unknown (Env samples)
# #   // // // Continuous (Ct) model compartment
# #   Y_qPCR_st ~ normal(mu_st,sigma_st);//Standards
# #   Y_qPCR_en ~ normal(mu_en,sigma_en);//Unknown (Env samples)
# #   // // // ddPCR
# #   W_st ~ binomial(U_st, inv_cloglog(omega_st)); //Standards
# #   W_en ~ binomial(U_en, inv_cloglog(omega_en)); //Unknown (Env samples)
# #   //
# #     // Priors
# #   // // qPCR
# #   // // // Bernoulli model
# #   phi_0 ~ normal(0, 2);
# #   phi_1 ~ normal(0, 2);
# #   // // // Continuous model
# #   beta_0 ~ normal(0, 3);
# #   beta_1 ~ normal(-3, 0.1);
# #   gamma_0 ~ normal(0, 0.1);
# #   gamma_1 ~ normal(-1, 0.1);
# #   C_qPCR ~ normal(0,3);
# #   // // ddPCR
# #   kappa_0 ~ normal(0, 1);
# #   kappa_1 ~ normal(0, 1);
# #   C_ddPCR ~ normal(0,3);
# # }
# 
# #  2.3 Run stan Model
# # Gledis Guri edited his 'Stan_qPCR_ddPCR_quantification_precision.stan' on the
# # 26-Feb-2025. The new version is placed in the 'stancode2' directory and is
# # named 'Stan_qPCR_ddPCR_quantification_precision_2.stan'.
# # He also provided this answer:
#  
# # Sorry for the late reply, too many things going on.
# # 
# # The R hat and ESS is indicating that the model didn't converge. You should not continue the next steps if the model didn't converge. That means that the parameters are not being estimated properly. It's very hard to judge what is going on without having the data but here's a couple of steps for troubleshooting:
# # 1. Try to increase the tree_depth by adding a line of code 'control = list(max_treedepth = 12),' after 'warmup = 2000,' condition when you run stanMod_2 line.
# # 2. Try to thin the data to just a couple of (good) samples and then incrementally increase the amount of samples. Try to include from the whole range of concentrations, especially where the probability of detection is about 50%. So to conclude, take 2 samples from 50% detection, 2 samples from <50% detection and 2 samples from >50% detection and run them as the unknown data (keep the standards the same as previously).
# # 3. Try to use only 1 assay at a time, this way you might identify which assay (or target) is outlier and doesn't follow the established relationship (the stan model).
# # 4. Check that the variability between plates is not high, so if you have the same target in multiple plates and if each plate is (very) different from one another this might lead the model having a hard time to converge. The way to check is to plot raw Ct values against known log(DNA conc). Or maybe just include samples from only 1 plate (in point 1).
# # 5. I made a new model (attached) that maybe can fix your problem, if your problem is in the probability of detection (Bernoulli) part of the model. I ran it and it works very smoothly, so give it a try, also using the steps 1 - 3. You can call it in stanMod_M2 line, just replace 'Stan_qPCR_ddPCR_quantification_precision.stan' with the new model name file.
# # 
# # Try these steps and you can share which parameters are still having trouble converging. You can run extract_param(stanMod_M2,c('phi_0','phi_1','beta_0','beta_1','gamma_0','gamma_1','kappa_0','kappa_1')). If you are using the new model (attached here) just remove 'phi_1' from the list of parameters.
# # 
# # Note that increasing the number of iterations won't help the model converge in this situation, so just keep the runs at around 5000 or slightly less.
# # 
# # Let me know how it goes. I have slightly more free time now (so will reply faster)
# #
# stanMod_M2 = stan(
#   #file = here::here('code','Stan_qPCR_ddPCR_quantification_precision.stan'),
#   file = here::here('stancode02','Stan_qPCR_ddPCR_quantification_precision_2.stan'),
#                   chains = 4,
#                   model_name='qPCR_vs_ddPCR_quantification_precision.stan',
#                   iter = 5000,
#                   warmup = 2000,
#                   control = list(max_treedepth = 12),
#                   # the 'stanMod_M2' run will not finish
#                   # with a 'stanMod_M2' that has data in it
#                   # without the initial values.
#                   # Without the in  initial values the model
#                   # is returned with an error saying that
#                   # 'Stan model 'qPCR_vs_ddPCR_quantification_precision.stan' does not contain samples.'
#                   init = initial_values,
#                   data = stan_data_M2)
# 
# #  2.4 Data for Figure 2
# 
# # Extract DNA concentration estimation from qPCR model
# q <- extract_param(stanMod_M2,"C_qPCR") 
# 
# # Extract DNA concentration estimation from ddPCR model
# d <- extract_param(stanMod_M2,"C_ddPCR") 
# 
# # Combine the data together
# diff <- 
#   d %>% 
#   dplyr::select(mean,`2.5%`,`97.5%`) %>% 
#   setNames(c("C_d","d2","d98")) %>%
#   cbind(.,q %>% 
#           dplyr::select(mean,`2.5%`,`97.5%`) %>% 
#           setNames(c("C_q","q2","q98"))) %>%
#   cbind(.,e_ddpcr %>% 
#           dplyr::filter(!duplicated(sp_st_idx)) %>% 
#           dplyr::arrange(sp_st_idx))
# 
# # 2.5 Plot Figure 2 - Quantification precision (qPCR vs ddPCR)
# 
# # store the diff data frame with the Quantification precision
# # in the outdir02 directory, with the species set number and 
# # concatenated species abbreviations appended to the file name
# write.csv(diff,
#           paste0(outdir02,"/diff_ddPCR_qPCR_quant_prec_",
#                  no_sbs,"_",
#                  sbs_spcs_cnc,
#                  ".csv"),
#           row.names = F)
# # source(here('Code','Surpressed','Figure_2.R'))
# # Figure_2
# 
# pp1 <-
#   diff %>%
#   ggplot()+
#   geom_errorbar(aes(y=10^C_d,
#                     x=10^C_q,
#                     ymin=10^d2,
#                     ymax=10^d98),
#                 color="grey")+
#   geom_errorbar(aes(y=10^C_d,
#                     x=10^C_q,
#                     xmin=10^q2,
#                     xmax=10^q98),
#                 color="grey")+
#   geom_point(aes(y=10^C_d,
#                  x=10^C_q,
#                  color=Species))+
#   scale_color_manual(values=
#                        # use the color range for the number
#                        # of species in the data frame
#                        # as inferred above
#                        c(clr04))+
#                        # c("tomato2", "deepskyblue2","orange2"))+
#   geom_abline(intercept = 0,slope=1,lty=2)+
#   theme_bw()+
#   scale_y_log10(labels=scientific_10,
#                 lim=c(1e-4,1e4),
#                 breaks=10^c(-4,-2,0,2,4))+
#   scale_x_log10(labels=scientific_10,
#                 lim=c(1e-5,2e2),
#                 breaks=10^c(-4,-2,0,2))+
#   facet_grid(~Species)+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
#         axis.title.x = element_text(size = 19),
#         axis.text.x = element_text(size=18),
#         axis.text.y=element_text(size=14))+
#   ylab("ddPCR modelled concentration \n(copies/µL)")+
#   xlab("qPCR modelled concentration (copies/µL)")
# 
# #pp1
# 
# # save the figure if the above 'bSaveFigures' is TRUE
# if(bSaveFigures==T){
#   ggsave(plot = pp1, 
#          # define the output filenmae by pasting together 
#          # qpcrrunno and qpcrrundate
#          filename = paste0(outdir02,
#         "/Fig02_01_ddPCR_and_qPCR_modelled_concentration",
#                            ".png"),
#          width=210*1.6,height=297*0.6,
#          #width=297,height=210,
#          units="mm",dpi=300)
# }
# 
# pp2 <-
#   diff %>% 
#   dplyr::filter(C_d>-2) %>%
#   dplyr::filter(C_q>-2) %>%
#   ggplot()+
#   geom_smooth(aes(y=10^(d98-d2),x=10^C_d,color=Species),
#               lty=2,se=T)+
#   geom_smooth(aes(y=10^(q98-q2),x=10^C_q,color=Species),
#               lty=1,se=T)+
#   scale_color_manual(values=
#                        # use the color range for the number
#                        # of species in the data frame
#                        # as inferred above
#                        c(clr04))+
#   # c("tomato2", "deepskyblue2","orange2"))+
#   scale_y_log10(labels=scientific_10,lim=c(1e0,1e3))+
#   scale_x_log10(labels=scientific_10,lim=c(1e-5,1e2),
#                 breaks=10^c(-4,-2,0,2))+
#   theme_bw()+
#   facet_grid(~Species)+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
#         axis.title.x = element_text(size = 19),
#         axis.text.x = element_text(size=18),
#         axis.text.y=element_text(size=14))+
#   ylab(expression("95% Credible Interval range"))+
#   xlab("Modelled concentration (copies/µL)")
# 
# 
# p1_leg <- ggplot() +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[1]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[2]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[3]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[4]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[5]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[6]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[7]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[8]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[9]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[10]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[11]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[12]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[13]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[14]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[15]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[16]), size = 1) +
#   geom_point(aes(x = NA, y = NA, 
#                  color = df_sp.sq$spcabr.nm[17]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[1]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[2]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[3]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[4]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[5]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[6]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[7]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[8]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[9]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[10]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[11]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[12]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[13]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[14]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[15]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[16]), size = 1) +
#   geom_line(aes(x = NA, y = NA, 
#                 color = df_sp.sq$spcabr.nm[17]), size = 1) +
#   theme_bw()+
#   scale_color_manual(
#     name='Assay',
#     # breaks=c("Gadus morhua",
#     #          "Clupea harengus",
#     #          "Pollachius virens"),
#     # values = c("tomato2", "deepskyblue2","orange"),
#     breaks=c(df_sp.sq$spcabr.nm),
#     values = c(clr04),
#     guide = guide_legend(
#       override.aes = list(
#         #lty = c(1,1,1),
#         #lty = c(1,1,1),
#         size = c( rep(3,Nspc) )))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16,face = "italic"))
# 
# p1_leg_leg <- cowplot::get_legend(p1_leg+
#                                     theme(legend.justification = 
#                                             c(0,0.1)))
# 
# p2_leg <- ggplot() +
#   geom_line(aes(x = NA, y = NA, color = "ddPCR"), size = 1) +
#   geom_line(aes(x = NA, y = NA, color = "qPCR"), size = 1) +
#   theme_bw()+
#   scale_color_manual(
#     name="Quantification precision \n(expressed as variance)",
#     breaks=c("ddPCR","qPCR"),
#     values = c("black","black"),
#     guide = guide_legend(
#       override.aes = list(lty = c(2,1),
#                           size = c(1,2)))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16,face = "italic"))
# p2_leg_leg <- cowplot::get_legend(p2_leg+
#                                     theme(legend.justification = 
#                                             c(0,1.0)))
# 
# legend <- cowplot::plot_grid(p1_leg_leg,p2_leg_leg,nrow = 2)
# p1 <- cowplot::plot_grid(pp1,pp2,nrow=2,
#                          align = "v",
#                          rel_heights = c(4,3),
#                          label_size = 20,
#                          labels = c("a","b"),
#                          label_x = -0.003)
# Figure_2 <- cowplot::plot_grid(p1,
#                                legend,
#                                ncol=2,
#                                rel_widths = c(4,1.0))
# Figure_2
# 
# # save the figure if the above 'bSaveFigures' is TRUE
# if(bSaveFigures==T){
#   ggsave(plot = Figure_2, 
#          # define the output filenmae by pasting together 
#          # qpcrrunno and qpcrrundate
#          filename = paste0(outdir02,
#       "/Fig02_02_ddPCR_and_qPCR_modelled_concentration",
#                            ".png"),
#          width=210*2.4,height=297*0.8,
#          #width=297,height=210,
#          units="mm",dpi=300)
# }
# 
# #_______________________________________________________________________________
# 
# 
# # 3 Two-step model
# # 3.1 Create stan data
# 
# # # Remove an outlier that has an influence on posteriors
# # e_qpcr_cm <- e_qpcr_cm %>% 
# #   dplyr::filter(!(Well=="A4"&
# #                     Sample_name=="2021624_10"&
# #                     Species =="Herring"&
# #                     Ct==42.56510544&
# #                     Conc==0.139179796&
# #                     plate=="Plate_G"))
# 
# stan_data_M3 <- list(
#   N_st_q = nrow(st_qpcr),
#   N_en_q = nrow(e_qpcr),
#   N_st_qp = nrow(st_qpcr_cm),
#   N_en_qp = nrow(e_qpcr_cm),
#   #
#   N_st_d = nrow(st_ddpcr),
#   N_en_d = nrow(e_ddpcr),
#   #
#   N_i = length(unique(st_qpcr$Species)),
#   N_ij = length(unique(e_ddpcr$sp_st_idx)),
#   #
#   i_qst_idx = st_qpcr$species_idx,
#   i_qen_idx = e_qpcr$species_idx,
#   ij_qen_idx = e_qpcr$sp_st_idx,
#   #
#   i_qst_p_idx = st_qpcr_cm$species_idx,
#   i_qen_p_idx = e_qpcr_cm$species_idx,
#   ij_qen_p_idx = e_qpcr_cm$sp_st_idx,
#   #
#   Z_qPCR_st = st_qpcr$pres,
#   Z_qPCR_en = e_qpcr$pres,
#   #C_qPCR_st = log10(st_qpcr$int_concentation),
#   C_qPCR_st = log10(st_qpcr$init_concentration),
#   #
#   Y_qPCR_st = st_qpcr_cm$Ct,
#   Y_qPCR_en = e_qpcr_cm$Ct,
#   #C_qPCR_st_cont = log10(st_qpcr_cm$int_concentation)
#   C_qPCR_st_cont = log10(st_qpcr_cm$init_concentration)
# );str(stan_data_M3)
# 
# # List of 20
# #  $ N_st_q        : int 277
# #  $ N_en_q        : int 1379
# #  $ N_st_qp       : int 242
# #  $ N_en_qp       : int 662
# #  $ N_st_d        : int 247
# #  $ N_en_d        : int 1052
# #  $ N_i           : int 3
# #  $ N_ij          : int 343
# #  $ i_qst_idx     : int [1:277] 1 1 1 1 1 1 1 1 1 1 ...
# #  $ i_qen_idx     : int [1:1379] 1 1 1 1 2 2 2 2 1 1 ...
# #  $ ij_qen_idx    : int [1:1379] 1 1 1 1 172 172 172 172 2 2 ...
# #  $ i_qst_p_idx   : int [1:242] 1 1 1 1 1 1 1 1 1 1 ...
# #  $ i_qen_p_idx   : int [1:662] 1 1 1 2 2 1 2 2 2 2 ...
# #  $ ij_qen_p_idx  : int [1:662] 1 1 2 173 173 3 174 175 175 175 ...
# #  $ Z_qPCR_st     : int [1:277] 1 1 1 1 1 1 1 1 1 1 ...
# #  $ Z_qPCR_en     : int [1:1379] 1 1 0 0 0 0 0 0 0 0 ...
# #  $ C_qPCR_st     : num [1:277] 6 5 4 3 2 1 0 6 4 3 ...
# #  $ Y_qPCR_st     : num [1:242] 17.2 20.3 24.3 26.9 30.7 ...
# #  $ Y_qPCR_en     : num [1:662] 43.9 46.2 49.6 40.6 40.4 ...
# #  $ C_qPCR_st_cont: num [1:242] 6 5 4 3 2 1 0 6 4 3 ...
# # 
# # 3.2 The math
# # 
# # 3.3 The stan code
# # 
# # data { 
# #   //Numbers of dimentions
# #   // // // qPCR
# #   int N_st_q; // Total number of observation in qPCR standard samples
# #   int N_en_q; // Total number of observation in qPCR environmental samples
# #   int N_st_qp; // Total number of observation in qPCR standard samples for only detected samples
# #   int N_en_qp; // Total number of observation in qPCR environmental samples for only detected samples
# #   // // // ddPCR
# #   int N_en_d; // Total number of observation in qPCR environmental samples
# #   // // Joined
# #   int N_i; // Number of species in both data
# #   int N_ij; // Number of species and stations in both data
# #   // 
# #     //Indexes
# #   // // qPCR
# #   // // // Binomial model
# #   array[N_st_q] int i_qst_idx; // Species index for qPCR standard samples
# #   array[N_en_q] int i_qen_idx; // Species index for qPCR environmental samples
# #   array[N_en_q] int ij_qen_idx; // Species and standard index for qPCR environmental samples
# #   // // // Continious model
# #   array[N_st_qp] int i_qst_p_idx; // Species index for qPCR standard samples
# #   array[N_en_qp] int i_qen_p_idx; // Species index for qPCR environmental samples
# #   array[N_en_qp] int ij_qen_p_idx; // Species and standard index for qPCR environmental samples
# #   // 
# #     // Data
# #   // // qPCR
# #   // // // Binomial model
# #   array[N_st_q] int Z_qPCR_st; // Presence/Absence response of qPCR standard data
# #   array[N_en_q] int Z_qPCR_en; // Presence/Absence response of qPCR environmental data
# #   array[N_st_q] real C_qPCR_st; // Known concentration (log10) in qPCR data
# #   // // // Continious model
# #   array[N_st_qp] real Y_qPCR_st; // Ct values of qPCR standard data for only detected samples
# #   array[N_en_qp] real Y_qPCR_en; // Ct values of qPCR environmental data for only detected samples
# #   array[N_st_qp] real C_qPCR_st_cont; // Known concentration (log10) in qPCR data for only detected samples
# # }
# # parameters {
# #   // Parameters
# #   // // qPCR
# #   // // // Bernoulli model
# #   vector[N_i] phi_0;
# #   vector[N_i] phi_1;
# #   // // // Continous model (joint)
# #   vector[N_i] beta_joint_0;
# #   vector[N_i] beta_joint_1;
# #   vector[N_i] gamma_joint_0;
# #   vector<upper=0>[N_i] gamma_joint_1;
# #   vector[N_ij] C_qPCR_joint; 
# #   // // // Continous model (single)
# #   vector[N_i] beta_cont_0;
# #   vector[N_i] beta_cont_1;
# #   vector[N_i] gamma_cont_0;
# #   vector<upper=0>[N_i] gamma_cont_1;
# #   vector[N_ij] C_qPCR_continuous;
# # }
# # transformed parameters{
# #   // Parameters
# #   // // qPCR
# #   // // // Bernoulli model
# #   vector[N_st_q] theta_st;
# #   vector[N_en_q] theta_un;
# #   // // // Continious model (joint)
# #   vector[N_st_qp] mu_st;
# #   vector[N_en_qp] mu_en;
# #   vector[N_st_qp] sigma_st;
# #   vector[N_en_qp] sigma_en;
# #   // // // Continious model (single)
# #   vector[N_st_qp] mu_st_nor;
# #   vector[N_en_qp] mu_en_nor;
# #   vector[N_st_qp] sigma_st_nor;
# #   vector[N_en_qp] sigma_en_nor;
# #   // 
# #     // Model TP
# #   // // qPCR model
# #   // // // Bernuli module model compartment
# #   // // // // // Standard
# #   for (i in 1:N_st_q){
# #     theta_st[i] = phi_0[i_qst_idx[i]] + (phi_1[i_qst_idx[i]] * C_qPCR_st[i]);
# #   }
# #   // // // // // Unknown
# #   for (i in 1:N_en_q){
# #     theta_un[i] = phi_0[i_qen_idx[i]] + (phi_1[i_qen_idx[i]] * C_qPCR_joint[ij_qen_idx[i]]);
# #   }
# #   // // // Continious model compartment (joint)
# #   // // // // Standard
# #   for (i in 1:N_st_qp){
# #     mu_st[i] = beta_joint_0[i_qst_p_idx[i]] + (beta_joint_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]);
# #     sigma_st[i] = exp(gamma_joint_0[i_qst_p_idx[i]]+(gamma_joint_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]));
# #   }
# #   // // // // Unknown
# #   for (i in 1:N_en_qp){
# #     mu_en[i] = beta_joint_0[i_qen_p_idx[i]] + (beta_joint_1[i_qen_p_idx[i]] * C_qPCR_joint[ij_qen_p_idx[i]]);
# #     sigma_en[i] = exp(gamma_joint_0[i_qen_p_idx[i]]+(gamma_joint_1[i_qen_p_idx[i]] * C_qPCR_joint[ij_qen_p_idx[i]]));
# #   }
# #   // // // Continious model compartment (single)
# #   // // // // Standard
# #   for (i in 1:N_st_qp){
# #     mu_st_nor[i] = beta_cont_0[i_qst_p_idx[i]] + (beta_cont_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]);
# #     sigma_st_nor[i] = exp(gamma_cont_0[i_qst_p_idx[i]]+(gamma_cont_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]));
# #   }
# #   // // // // Unknown
# #   for (i in 1:N_en_qp){
# #     mu_en_nor[i] = beta_cont_0[i_qen_p_idx[i]] + (beta_cont_1[i_qen_p_idx[i]] * C_qPCR_continuous[ij_qen_p_idx[i]]);
# #     sigma_en_nor[i] = exp(gamma_cont_0[i_qen_p_idx[i]]+(gamma_cont_1[i_qen_p_idx[i]] * C_qPCR_continuous[ij_qen_p_idx[i]]));
# #   }
# # }
# # model {
# #   // Model 
# #   // // qPCR
# #   // // // Bernoulli model
# #   Z_qPCR_st ~ bernoulli(inv_logit(theta_st)); //Standards
# #   Z_qPCR_en ~ bernoulli(inv_logit (theta_un)); //Unknown (Env samples)
# #   // // // Continuous (Ct) model compartment (joint)
# #   Y_qPCR_st ~ normal(mu_st,sigma_st);//Standards
# #   Y_qPCR_en ~ normal(mu_en,sigma_en);//Unknown (Env samples)
# #   // // // Continuous (Ct) model compartment (single)
# #   Y_qPCR_st ~ normal(mu_st_nor,sigma_st_nor);//Standards
# #   Y_qPCR_en ~ normal(mu_en_nor,sigma_en_nor);//Unknown (Env samples)
# #   // 
# #     // Priors
# #   // // qPCR
# #   // // // Bernoulli model
# #   phi_0 ~ normal(0, 2);
# #   phi_1 ~ normal(0, 2);
# #   // // // Continious model (joint)
# #   beta_joint_0 ~ normal(0, 3);
# #   beta_joint_1 ~ normal(-3, 0.1);
# #   gamma_joint_0 ~ normal(0, 0.1);
# #   gamma_joint_1 ~ normal(-1, 0.1);
# #   C_qPCR_joint ~ normal(0,1);
# #   // // // Continious model (single)
# #   beta_cont_0 ~ normal(0, 3);
# #   beta_cont_1 ~ normal(-3, 0.1);
# #   gamma_cont_0 ~ normal(0, 0.1);
# #   gamma_cont_1 ~ normal(-1, 0.1);
# #   C_qPCR_continuous ~ normal(0,1);
# # }
# 
# # 3.4 Run stan Model
# 
# stanMod_M3 = stan(file = here::here('code','Stan_qPCR_two_step_model.stan'),
#                   chains = 4,
#                   model_name='qPCR_two_step_model',
#                   iter = 5000,
#                   warmup=2000,
#                   data = stan_data_M3)
# 
# # Warning: Tail Effective Samples Size (ESS) is too low, 
# # indicating posterior variances and tail quantiles may be unreliable.
# # Running the chains for more iterations may help. See
# # https://mc-stan.org/misc/warnings.html#tail-ess
# # 
# # 3.5 Data for Figure 3
# 
# q_bin <-
#   extract_param(stanMod_M3,"C_qPCR_joint") %>% 
#   cbind(.,e_qpcr %>% 
#           dplyr::filter(!duplicated(sp_st_idx)) %>% 
#           dplyr::arrange(sp_st_idx) %>% 
#           dplyr::select(-Ct,-Conc,-pres)) %>% 
#   dplyr::select(mean,`97.5%`,`2.5%`,Species) %>% 
#   setNames(c("C_qPCR_joint","q98","q2","Species")) %>% 
#   dplyr::slice(e_qpcr_cm %>%
#                  dplyr::filter(!duplicated(sp_st_idx)) %>%
#                  dplyr::pull(sp_st_idx) %>%
#                  sort())
# q_nor <-
#   extract_param(stanMod_M3,"C_qPCR_continuous") %>% 
#   dplyr::slice(e_qpcr_cm %>% 
#                  dplyr::filter(!duplicated(sp_st_idx)) %>% 
#                  dplyr::pull(sp_st_idx) %>% 
#                  sort()) %>% 
#   cbind(.,e_qpcr_cm %>% 
#           dplyr::filter(!duplicated(sp_st_idx)) %>% 
#           dplyr::arrange(sp_st_idx) %>% 
#           dplyr::select(-Ct,-Conc,-pres)) %>% 
#   dplyr::select(mean,`97.5%`,`2.5%`,Species) %>% 
#   setNames(c("C_qPCR_continuous","q98","q2","Species"))
# 
# # the 'diff' data frame from the webpage
# diff <- cbind(q_bin %>% 
#                 setNames(c("C_qPCR_joint",
#                            "q98_bin",
#                            "q2_bin",
#                            "Species")) %>% 
#                 dplyr::select(-Species),
#               q_nor %>% setNames(c("C_qPCR_continuous",
#                                    "q98_nor",
#                                    "q2_nor",
#                                    "Species"))) %>% 
#   dplyr::mutate(diff=(q98_nor-q2_nor)-(q98_bin-q2_bin))
# 
# # Other version of 'diff' that creates the 'C_bin' column
# diff <- cbind(q_bin %>% setNames(c("C_bin",
#                                    "q98_bin",
#                                    "q2_bin",
#                                    "Species")) %>% 
#                 dplyr::select(-Species),
#               q_nor %>% setNames(c("C_nor",
#                                    "q98_nor",
#                                    "q2_nor",
#                                    "Species"))) %>% 
#   dplyr::mutate(diff=(q98_nor-q2_nor)-(q98_bin-q2_bin))
# 
# # store the data frame with the Two-step qPCR model 
# # quantification precision
# # inferred from the stan model 3. Store it as a data frame in the 
# # outdir02 directory, with the species number
# # and the concatenated species nmawes in the file name
# write.csv(diff,
#           paste0(outdir02,"/diff_Two_step_qPCR_mod_quant_",
#                  no_sbs,"_",
#                  sbs_spcs_cnc,".csv"),
#           row.names=F)
# # 3.6 Plot Figure 3 - Two-step qPCR model quantification precision
# 
# pp1 <-
#   diff %>%
#   ggplot()+
#   geom_point(aes(y=I(q98_bin-q2_bin),
#                  x=10^C_bin,color=Species),
#              size=2,shape=16)+
#   geom_point(aes(y=I(q98_nor-q2_nor),
#                  x=10^C_nor,color=Species),
#              size=2,shape=3)+
#   scale_color_manual(values=
#                        # use the color range inferred above
#                        c(clr04))+
#                        #c("tomato2", "deepskyblue2","orange2"))+
#   scale_x_log10()+
#   theme_bw()+
#   facet_grid(~Species)+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         plot.margin = margin(0.1, 0.1, 0.5, 1, "cm"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y=element_text(size=14))+
#   ylab(expression(atop("Quantification variance \n(measured as 95% CI)",
#                        paste("log"[10]*"(copies/µL)"))))
# pp1 
# 
# pp2 <-
#   diff %>%
#   ggplot()+
#   geom_smooth(aes(y=diff,
#                   x=10^C_bin,color=Species),
#               lty=2,se=F)+
#   geom_point(aes(y=diff,
#                  x=10^C_bin,color=Species),
#              size=2,shape=4)+
#   scale_color_manual(values=
#                        # use the color range inferred above
#                        c(clr04))+
#   #c("tomato2", "deepskyblue2","orange2"))+
#   theme_bw()+
#   scale_x_log10(labels=scientific_10)+
#   facet_grid(~Species)+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         plot.margin = margin(0.1, 0.1, 0.5, 1, "cm"),
#         axis.title.x = element_text(size = 19),
#         axis.text.x = element_text(size=18),
#         axis.text.y=element_text(size=14))+
#   ylab(expression(atop(
#     "Quantification variance difference \n     (measured as CI difference)",
#     paste("(CI"[continuous]*" - CI"[logistic*" + "* continuous]*")"))))+
#   ylim(-0.3,0.8)+
#   xlab("Modelled qPCR concentration (copies/µL)")
# 
# pp2
# 
# 
# p1 <- cowplot::plot_grid(pp1,pp2,
#                          nrow = 2,
#                          align = "v",
#                          rel_heights = c(3.5,4),
#                          labels = c("a","b"),
#                          label_size = 20, label_x = 0.08)
# 
# p1_leg <- ggplot() +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[1]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[2]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[3]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[4]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[5]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[6]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[7]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[8]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[9]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[10]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[11]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[12]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[13]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[14]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[15]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[16]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[17]), size = 1) +
#   
#   theme_bw()+
#   scale_color_manual(
#     name='Species',
#     breaks=
#       # c("Gadus morhua", 
#       #   "Clupea harengus",
#       #   "Pollachius virens"),
#       c(df_sp.sq$spcabr.nm),
#     values = 
#       # use the color range inferred above
#       c(clr04))+
#   #c("tomato2", "deepskyblue2","orange2"))+
#     # guide = guide_legend(
#     #   override.aes = list(
#     #     #lty = c(1,1,1),
#     #     size = c( rep(2,Nspc) ))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16,face = "italic"))
# 
# p2_leg <- ggplot() +
#   geom_point(aes(x = NA, y = NA, color = "Joint model"), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = "Continuous model"), size = 1) +
#   theme_bw()+
#   scale_color_manual(
#     name='Model used',
#     breaks=c("Joint model", "Continuous model"),
#     values = c("black", "black"),
#     guide = guide_legend(
#       override.aes = list(shape = c(16,3),
#                           size = c(2,2)))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16))
# 
# p3_leg <- ggplot() +
#   geom_point(aes(x = NA, 
#                  y = NA, 
#                  color = "Precision difference"), 
#              size = 1) +
#   theme_bw()+
#   scale_color_manual(
#     name="",
#     breaks=c("Precision difference"),
#     values = c("black"),
#     guide = guide_legend(
#       override.aes = list(shape = c(4),
#                           size = c(2)))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=1),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16))
# 
# 
# p1_leg_leg <- cowplot::get_legend(p1_leg+
#                                     theme(legend.justification = 
#                                             c(0,-2.8)))
# p2_leg_leg <- cowplot::get_legend(p2_leg+
#                                     theme(legend.justification = 
#                                             c(0,-0.4)))
# p3_leg_leg <- cowplot::get_legend(p3_leg+
#                                     theme(legend.justification = 
#                                             c(0,2.2)))
# 
# legend <-
#   cowplot::plot_grid(p1_leg_leg,
#                      p2_leg_leg,
#                      p3_leg_leg,
#                      nrow = 4)
# Figure_3 <- cowplot::plot_grid(p1,
#                                legend,
#                                ncol=2,
#                                rel_widths = c(4,2.8))
# Figure_3
# 
# # save the figure if the above 'bSaveFigures' is TRUE
# if(bSaveFigures==T){
#   ggsave(plot = Figure_3, 
#          # define the output filenmae by pasting together 
#          # qpcrrunno and qpcrrundate
#          filename = paste0(outdir02,
#       "/Fig03_01_ddPCR_and_qPCR_Quantification_variance_difference_",
#       no_sbs,"_",sbs_spcs_cnc,
#                            ".png"),
#          width=210*1.6,height=297*0.8,
#          #width=297,height=210,
#          units="mm",dpi=300)
# }
# #source(here('Code','Surpressed','Figure_3.R'))
# #Figure_3
# 
# # 4 Lower threshold
# # 4.1 The math
# 
# ## Extract the data from stan Model 2
# 
# kappa_0 <- extract_param(stanMod_M2,"kappa_0") %>%
#   cbind(.,e_ddpcr %>%
#           dplyr::filter(!duplicated(species_idx)) %>%
#           dplyr::arrange(species_idx) %>%
#           dplyr::select(Species,species_idx))
# kappa_1 <- extract_param(stanMod_M2,"kappa_1") %>%
#   cbind(.,e_ddpcr %>%
#           dplyr::filter(!duplicated(species_idx)) %>%
#           dplyr::arrange(species_idx) %>%
#           dplyr::select(Species,species_idx))
# 
# nrep_st <- 1
# nrep_end <- 10
# repl <- rep(c(nrep_st:nrep_end),3)
# 
# minimum_threshold <- repl %>% as.data.frame() %>%
#   cbind(.,kappa_0 %>% 
#           dplyr::pull(mean) %>%
#           rep(.,each=nrep_end)) %>%
#   cbind(.,kappa_1 %>% 
#           dplyr::pull(mean) %>% 
#           rep(.,each=nrep_end)) %>%
#   cbind(.,kappa_0 %>%
#           dplyr::pull(Species) %>% 
#           rep(.,each=nrep_end)) %>%
#   setNames(c("replicates","kappa_0","kappa_1","Species")) %>%
#   dplyr::mutate(C_lt=(cloglog(1/(20000*replicates))-
#                         kappa_0)/kappa_1) %>%
#   dplyr::mutate(C_ut=(cloglog(((20000*replicates)-1)/
#                                 (20000*replicates))-kappa_0)/kappa_1)
# # other minimum_threshold data frame from R script
# minimum_threshold <- repl %>% 
#   as.data.frame() %>%
#   cbind(.,kappa_0 %>% 
#           dplyr::pull(mean) %>% 
#           rep(.,each=nrep_end)) %>%
#   cbind(.,kappa_1 %>% 
#           dplyr::pull(mean) %>% 
#           rep(.,each=nrep_end)) %>%
#   cbind(.,kappa_0 %>% 
#           dplyr::pull(Species) %>% 
#           rep(.,each=nrep_end)) %>%
#   setNames(c("replicates","intercept","slope","Species")) %>%
#   dplyr::mutate(C_th=(cloglog(1/(20000*replicates))-
#                         intercept)/slope) %>%
#   dplyr::mutate(C_th_up=(cloglog(((20000*replicates)-1)/
#                                    (20000*replicates))-intercept)/slope)
# 
# # store the data frame with the minimum threshold
# # in the outdir02 directory, with the species number
# # and the concatenated species names in the file name
# write.csv(minimum_threshold,
#           paste0(outdir02,"/df_minimum_threshold_",
#                  no_sbs,"_",
#                  sbs_spcs_cnc,".csv"),
#           row.names=F)
# # 4.2 Plot Figure 4 - ddPCR minimum threshold
# 
# pp1 <-
#   minimum_threshold %>%
#   ggplot()+
#   geom_line(aes(y=C_th,x=replicates,color=Species))+
#   geom_point(aes(y=C_th,x=replicates,color=Species))+
#   scale_x_continuous(breaks=c(1,3,5,7,10))+
#   scale_y_continuous(breaks=seq(-2.5,-1,by=0.2))+
#   scale_color_manual(
#     name='Species',
#     breaks=
#       # c("Gadus morhua", 
#       #   "Clupea harengus",
#       #   "Pollachius virens"),
#       c(df_sp.sq$spcabr.nm),
#     values = 
#       # use the color range inferred above
#       c(clr04))+
#   theme_bw()+
#   facet_grid(~Species)+
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(size = 19),
#         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
#         axis.title.x = element_text(size = 19),
#         axis.text.x = element_text(size=14),
#         axis.text.y=element_text(size=14))+
#   ylab(expression(atop("Lower threshold",
#                        paste("Log"[10]*" concentration (copies/µL)"))))+
#   xlab("Number of replicates")
# 
# p1_leg <- ggplot() +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[1]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[2]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[3]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[4]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[5]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[6]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[7]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[8]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[9]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[10]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[11]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[12]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[13]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[14]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[15]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[16]), size = 1) +
#   geom_point(aes(x = NA, y = NA, color = df_sp.sq$spcabr.nm[17]), size = 1) +
#   
#   
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[1]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[2]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[3]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[4]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[5]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[6]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[7]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[8]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[9]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[10]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[11]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[12]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[13]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[14]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[15]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[16]), size = 1) +
#   geom_line(aes(x = NA, y = NA,color = df_sp.sq$spcabr.nm[17]), size = 1) +
#   
#   theme_bw()+
#   scale_color_manual(
#     name='Species',
#     breaks=
#       # c("Gadus morhua", 
#       #   "Clupea harengus",
#       #   "Pollachius virens"),
#       c(df_sp.sq$spcabr.nm),
#     values = 
#       # use the color range inferred above
#       c(clr04))+
#   #c("tomato2", "deepskyblue2","orange2"))+
#     # guide = guide_legend(
#     #   override.aes = list(lty = c(rep(1,Nspc)),
#     #                       size = c(rep(1,Nspc)))) +
#   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
#         legend.key.width = unit(1.2,"cm"),
#         legend.title = element_text(size=18),
#         legend.key.height = unit(0.8, 'cm'),
#         legend.text = element_text(size=16,face = "italic"))
# 
# p1_leg_leg <- cowplot::get_legend(p1_leg+
#                                     theme(legend.justification = 
#                                             c(0,0.5)))
# 
# 
# Figure_4 <- cowplot::plot_grid(pp1,
#                                p1_leg_leg,
#                                nrow=1,
#                                ncol=2,
#                                rel_widths = c(4,1))
# Figure_4
# 
# #source(here('code','Surpressed','Figure_4.R'))
# # save the figure if the above 'bSaveFigures' is TRUE
# if(bSaveFigures==T){
#   ggsave(plot = Figure_4, 
#          # define the output filenmae by pasting together 
#          # qpcrrunno and qpcrrundate
#          filename = paste0(outdir02,"/Fig04_01_ddPCR_and_qPCR_Lower_threshold_",
#                            no_sbs,"_",sbs_spcs_cnc,
#                            ".png"),
#          width=210*1.6,height=297*0.8,
#          #width=297,height=210,
#          units="mm",dpi=300)
# }
# #
# 
# # end iteration over the sets of species
# }
# 
# #______________________________________________________________________
# #______________________________________________________________________
# #______________________________________________________________________
# # 
# # # 5 Supplementary plots
# # # 
# # # Load additional data
# # 
# # #qpcr_optimisation_cod_herring_plate <- read.csv(here('Data','qPCR_opt_Cod_Herring.csv'))
# # qpcr_optimisation_cod_herring_plate <- qPCR_opt_Cod_Herring.csv
# # #qpcr_optimisation_cod_saithe_plate <- read.csv(here('Data','qPCR_opt_Cod_Saithe.csv'))
# # qpcr_optimisation_cod_saithe_plate <- qPCR_opt_Cod_Saithe.csv
# # 
# # #ddpcr_amplitude_cod_herring_standard_samples <- read.csv(here('Data','ddpcr_amplitude_cod_herring_standard_samples.csv'))
# # ddpcr_amplitude_cod_herring_standard_samples <- ddpcr_amplitude_cod_herring_standard_samples.csv
# # #ddpcr_amplitude_cod_saithe_standard_samples <- read.csv(here('Data','ddpcr_amplitude_cod_saithe_standard_samples.csv'))
# # ddpcr_amplitude_cod_saithe_standard_samples <- ddpcr_amplitude_cod_saithe_standard_samples.csv
# # # 5.1 Figure S1 - Optimization of qPCR plates
# # 
# # opt_dat <- rbind(qpcr_optimisation_cod_saithe_plate,
# #                  qpcr_optimisation_cod_herring_plate) %>% 
# #   dplyr::filter(SuperMix!="OUT") %>% 
# #   dplyr::mutate(Conc=if_else(Sample=="ST1",6,0)) %>% 
# #   dplyr::mutate(Conc=if_else(Sample=="ST3",4,Conc)) %>% 
# #   dplyr::mutate(Conc=if_else(Sample=="ST5",2,Conc)) %>% 
# #   dplyr::mutate(Conc=if_else(Sample=="Blank",-Inf,Conc)) %>% 
# #   dplyr::mutate(SM=substr(SuperMix,start = 3,stop = 3))
# # 
# # pp1 <-
# #   opt_dat %>% 
# #   dplyr::filter(Conc!=-Inf) %>% 
# #   ggplot(aes(x=SM,y=Ct,color=Species))+
# #   geom_boxplot()+
# #   theme_bw()+
# #   facet_wrap(~Conc ~ Species,scales="free_y",ncol=3,
# #              labeller=label_bquote(cols = 10 ^ .(Conc) ~ "(copies/µL)"))+
# #   scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
# #   ylab("Cycle threshold (Ct)")+
# #   xlab("SuperMix (SM)")+
# #   theme(legend.position = "none",
# #         strip.background = element_blank(),
# #         strip.text.x = element_text(size=17),
# #         axis.title.y = element_text(size = 19),
# #         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
# #         axis.title.x = element_text(size = 19),
# #         axis.text.x = element_text(size=15),
# #         axis.text.y=element_text(size=15))
# # 
# # p1_leg <-
# #   opt_dat %>% dplyr::filter(SuperMix!="OUT"&Conc!=-Inf) %>% 
# #   ggplot(aes(x=SuperMix,y=Ct,color=Species))+
# #   geom_boxplot()+
# #   theme_bw()+
# #   scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
# #   theme(legend.title = element_text(size = 20),
# #         legend.text = element_text(size = 16))  # Adjust the size as needed
# # 
# # mat <- matrix(paste0("SM", c(1:3,5:7,9,8,4)), 
# #               nrow = 3, byrow = TRUE)
# # 
# # df <- expand.grid(row = 1:3, col = 1:3)
# # df$label <- as.vector(mat)
# # df$F_p <- rep(c(2.5,10,17),3)
# # df$R_p <- rep(c(2.5,10,17),each=3)
# # 
# # pp2 <-
# #   ggplot(df, aes(x = col, y = row, label = label, fill = label)) +
# #   geom_tile(color = "white", size = 1) +
# #   geom_text(size = 5, color = "black") +
# #   geom_text(aes(x = 0.2, y = row, label = F_p), 
# #             hjust = 0.1,size=5) +
# #   geom_text(aes(x = col, y = 0.15, label = R_p), 
# #             vjust = -0.5,size=5) +
# #   scale_fill_manual(values = rep("grey", 9)) +
# #   ggtitle("SuperMix (SM) \nconcentration") +  # Add your desired title here
# #   theme_minimal() +
# #   ylab("Forward primer \nconcentration (nM)")+
# #   xlab("Reverse primer \nconcentration (nM)")+
# #   theme(
# #     legend.position = "none",
# #     axis.title.y = element_text(size = 14),
# #     axis.title.x = element_text(size = 14),
# #     axis.title = element_blank(),
# #     plot.margin = margin(1, 0, 0.1, 0.5, "cm"),
# #     axis.text = element_blank(),
# #     axis.ticks = element_blank(),
# #     panel.grid = element_blank(),
# #     plot.title = element_text(hjust = 0.5,size = 20)
# #   )
# # 
# # p1_leg_leg <- cowplot::get_legend(p1_leg+
# #                                     theme(legend.justification = 
# #                                             c(0.4,0.0)))
# # 
# # p1 <-
# #   cowplot::plot_grid(p1_leg_leg,pp2, nrow = 3)
# # 
# # Figure_S1 <- cowplot::plot_grid(pp1,p1,rel_widths = c(5,1.5))
# # 
# # Figure_S1
# # #source(here('code','Surpressed','Figure_4.R'))
# # # save the figure if the above 'bSaveFigures' is TRUE
# # if(bSaveFigures==T){
# #   ggsave(plot = Figure_S1, 
# #          # define the output filenmae by pasting together 
# #          # qpcrrunno and qpcrrundate
# #          filename = paste0(outdir02,"/FigS01_ddPCR_and_qPCR_cycle_threshold",
# #                            ".png"),
# #          width=210*1.6,height=297*0.8,
# #          #width=297,height=210,
# #          units="mm",dpi=300)
# # }
# # 
# # # 5.2 Figure S2 - ddPCR droplet proportion
# # 
# # pp1 <-
# #   ggplot(data=st_ddpcr %>%
# #            group_by(Sample_name, Species) %>%
# #            summarise(nominal_c = mean(int_concentation),
# #                      pos_agg = sum(Positives),
# #                      tot_agg = sum(Tot_drop),
# #                      sp_idx = mean(species_idx),
# #                      Conc = mean(Conc))) +
# #   geom_point(aes(x=Conc, y=cloglog(pos_agg/tot_agg),
# #                  color=Species),size=2)+
# #   geom_point(aes(x=nominal_c, y=cloglog(pos_agg/tot_agg),
# #                  color=Species),shape=23,size=2)+
# #   geom_smooth(aes(x=Conc, y=cloglog(pos_agg/tot_agg),
# #                   color=Species),
# #               method=lm, se=F, fullrange=T,lty=3,size=0.5)+
# #   stat_smooth(aes(x=nominal_c, y=cloglog(pos_agg/tot_agg),
# #                   color=Species),
# #               method = lm, formula = y ~ poly(x, 3),
# #               lty=1,size=0.5, se=F)+
# #   theme_bw()+
# #   scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
# #   facet_grid(~Species)+
# #   scale_x_log10(labels=scientific_10,
# #                 lim=c(1e-3,1e5),
# #                 breaks=c(10^seq(-3,5,by=2)))+
# #   theme(legend.position = "none",
# #         strip.background = element_blank(),
# #         strip.text.x = element_blank(),
# #         axis.title.y = element_text(size = 19),
# #         plot.margin = margin(0.1, 0, 0.3, 0, "cm"),
# #         axis.title.x = element_text(size = 19),
# #         axis.text.x = element_text(size = 15),
# #         axis.text.y=element_text(size=15))+
# #   ylab(expression(atop("cloglog("*omega*")",
# #                        paste(omega*" = positive droplets / total droplets"))))+
# #   xlab("Nominal concentration (copies/µL)")
# # 
# # p1_leg <- ggplot() +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Gadus morhua"), size = 1) +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Clupea harengus"), size = 1) +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Pollachius virens"), size = 1) +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "Gadus morhua"), size = 1) +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "Clupea harengus"), size = 1) +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "Pollachius virens"), size = 1) +
# #   theme_bw()+
# #   scale_color_manual(
# #     name='Assay',
# #     breaks=c("Gadus morhua", "Clupea harengus","Pollachius virens"),
# #     values = c("tomato2", "deepskyblue2","orange"),
# #     guide = guide_legend(
# #       override.aes = list(lty = c(1,1,1),
# #                           size = c(3,3,3)))) +
# #   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
# #         legend.key.width = unit(1.2,"cm"),
# #         legend.title = element_text(size=18),
# #         legend.key.height = unit(0.8, 'cm'),
# #         legend.text = element_text(size=16,face = "italic"))
# # 
# # p1_leg_leg <- cowplot::get_legend(p1_leg+
# #                                     theme(legend.justification = 
# #                                             c(0,0.0)))
# # 
# # p2_leg <- ggplot() +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "ddPCR (in-built estimation)"),size = 1) +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "Nominal concentration"), size = 1) +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "ddPCR (in-built estimation)"), size = 1,shape=23) +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Nominal concentration"), size = 1) +
# #   theme_bw()+
# #   scale_color_manual(
# #     name='Method of measurement',
# #     breaks=c("ddPCR (in-built estimation)","Nominal concentration"),
# #     values = c("black","black"),
# #     guide = guide_legend(
# #       override.aes = list(lty = c(3,1),
# #                           size = c(3,3),
# #                           shape= c(19,23)))) +
# #   theme(plot.margin = margin(1, 0.1, 0, 0, "cm"),
# #         legend.key.width = unit(1.2,"cm"),
# #         legend.title = element_text(size=18),
# #         legend.key.height = unit(0.8, 'cm'),
# #         legend.text = element_text(size=16))
# # p2_leg_leg <- cowplot::get_legend(p2_leg+
# #                                     theme(legend.justification = 
# #                                             c(0,1.0)))
# # 
# # legend <- cowplot::plot_grid(p1_leg_leg, 
# #                              p2_leg_leg, 
# #                              nrow = 2)
# # 
# # Figure_S2 <- cowplot::plot_grid(pp1,
# #                                 legend,
# #                                 nrow=1,
# #                                 rel_widths = c(4,1.1))
# # 
# # Figure_S2
# # 
# # # save the figure if the above 'bSaveFigures' is TRUE
# # if(bSaveFigures==T){
# #   ggsave(plot = Figure_S2, 
# #          # define the output filenmae by pasting together 
# #          # qpcrrunno and qpcrrundate
# #          filename = paste0(outdir02,"/FigS02_ddPCR_and_qPCR_Nominal_concentration",
# #                            ".png"),
# #          width=210*1.6,height=297*0.8,
# #          #width=297,height=210,
# #          units="mm",dpi=300)
# # }
# # # 
# # # 5.3 Figure S3 - Possion statistics between assays
# # 
# # kappa_0 <- extract_param(stanMod_M2,"kappa_0") %>%
# #   cbind(.,e_ddpcr %>%
# #           dplyr::filter(!duplicated(species_idx)) %>%
# #           dplyr::arrange(species_idx) %>%
# #           dplyr::select(Species,species_idx))
# # kappa_1 <- extract_param(stanMod_M2,"kappa_1") %>%
# #   cbind(.,e_ddpcr %>%
# #           dplyr::filter(!duplicated(species_idx)) %>%
# #           dplyr::arrange(species_idx) %>%
# #           dplyr::select(Species,species_idx))
# # 
# # conc_st <- -2
# # conc_end <- 5
# # repl <- rep(seq(conc_st,conc_end,by=0.1),3)
# # logist_diff <-
# #   repl %>% as.data.frame() %>%
# #   cbind(.,kappa_0 %>% 
# #           dplyr::pull(mean) %>% 
# #           rep(.,each=length(repl)/3)) %>%
# #   cbind(.,kappa_1 %>% 
# #           dplyr::pull(mean) %>% 
# #           rep(.,each=length(repl)/3)) %>%
# #   cbind(.,kappa_0 %>% 
# #           dplyr::pull(Species) %>% 
# #           rep(.,each=length(repl)/3)) %>%
# #   setNames(c("C","intercept","slope","Species")) %>%
# #   dplyr::mutate(omega=intercept+(slope*C))
# # 
# # sim <- seq(-2,4.9,by=0.1) %>% 
# #   as.data.frame() %>% 
# #   setNames("C") %>%
# #   dplyr::mutate(omega=inv.cloglog(-7.07+(2.3*C)))
# # 
# # pp1 <-
# #   logist_diff %>%
# #   ggplot()+
# #   geom_line(aes(y=inv.cloglog(omega),
# #                 x=10^C,color=Species),lwd=1)+
# #   geom_line(data=sim,aes(y=(omega),
# #                          x=10^C),color="black",lty=2,lwd=0.3)+
# #   scale_x_log10(labels=scientific_10,
# #                 breaks=c(10^seq(conc_st,conc_end)))+ #
# #   scale_y_sqrt()+
# #   scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
# #   theme_bw()+
# #   theme(legend.position = "none",
# #         strip.background = element_blank(),
# #         strip.text.x = element_blank(),
# #         axis.title.y = element_text(size = 19),
# #         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
# #         axis.title.x = element_text(size = 19),
# #         axis.text.x = element_text(size=14),
# #         axis.text.y=element_text(size=14))+
# #   ylab(expression("Proportion of positive droplets"))+
# #   xlab("Nominal concentration (copies/µL)")
# # 
# # 
# # p1_leg <- ggplot() +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "Gadus morhua"), size = 1) +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "Clupea harengus"), size = 1) +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "Pollachius virens"), size = 1) +
# #   theme_bw()+
# #   scale_color_manual(
# #     name='Bayesian estimated\nassay',
# #     breaks=c("Gadus morhua", "Clupea harengus","Pollachius virens"),
# #     values = c("tomato2", "deepskyblue2","orange"),
# #     guide = guide_legend(
# #       override.aes = list(lty = c(1,1,1),
# #                           size = c(2,2,2)))) +
# #   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
# #         legend.key.width = unit(1.2,"cm"),
# #         legend.title = element_text(size=18),
# #         legend.key.height = unit(0.8, 'cm'),
# #         legend.text = element_text(size=16,face = "italic"))
# # 
# # p1_leg_leg <- cowplot::get_legend(p1_leg+
# #                                     theme(legend.justification = 
# #                                             c(0,0.1)))
# # p2_leg <- ggplot() +
# #   geom_line(aes(x = NA, y = NA, 
# #                 color = "ddPCR (in-built estimation)"), 
# #             size = 1) +
# #   theme_bw()+
# #   scale_color_manual(
# #     name='Poisson statistics',
# #     breaks=c("ddPCR (in-built estimation)"),
# #     values = c("black"),
# #     guide = guide_legend(
# #       override.aes = list(lty = c(2),
# #                           size = c(1)))) +
# #   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
# #         legend.key.width = unit(1.2,"cm"),
# #         legend.title = element_text(size=18),
# #         legend.key.height = unit(0.8, 'cm'),
# #         legend.text = element_text(size=16,face = "italic"))
# # 
# # p2_leg_leg <- cowplot::get_legend(p2_leg+
# #                                     theme(legend.justification = 
# #                                             c(0,1.0)))
# # 
# # legend <- cowplot::plot_grid(p1_leg_leg,
# #                              p2_leg_leg,
# #                              nrow = 2)
# # Figure_S3 <- cowplot::plot_grid(pp1,
# #                                 legend,nrow=1,
# #                                 ncol=2,
# #                                 rel_widths = c(4,1.3))
# # 
# # Figure_S3
# # 
# # # save the figure if the above 'bSaveFigures' is TRUE
# # if(bSaveFigures==T){
# #   ggsave(plot = Figure_S3, 
# #          # define the output filenmae by pasting together 
# #          # qpcrrunno and qpcrrundate
# #          filename = paste0(outdir02,"/FigS03_Poisson_statistics",
# #                            ".png"),
# #          width=210*1.6,height=297*0.8,
# #          #width=297,height=210,
# #          units="mm",dpi=300)
# # }
# # # 
# # 
# # # 5.4 Figure S4 - Assay performance amplitude
# # 
# # amp_dat <- rbind(ddpcr_amplitude_cod_herring_standard_samples,
# #                  ddpcr_amplitude_cod_saithe_standard_samples)
# # 
# # pp1 <-
# #   amp_dat %>% 
# #   ggplot() +
# #   geom_jitter(aes(y=MeanAmplitudeOfPositives,
# #                   x=Species,shape = assay),
# #               size=2,color="forestgreen")+
# #   geom_boxplot(aes(x = Species, 
# #                    y = MeanAmplitudeOfPositives),fill="forestgreen",
# #                alpha=0.5)+
# #   geom_jitter(aes(y=MeanAmplitudeOfNegatives,
# #                   x=Species,shape = assay),
# #               size=2,color="tomato2")+
# #   geom_boxplot(aes(x = Species, 
# #                    y = MeanAmplitudeOfNegatives),
# #                fill="tomato2",
# #                alpha=0.5)+
# #   theme_bw()+
# #   scale_shape_manual(values = c(17, 19)) +
# #   labs(x = "Species",
# #        y = "Mean amplitude of droplets")+
# #   theme(legend.position = "none",
# #         strip.background = element_blank(),
# #         strip.text.x = element_blank(),
# #         axis.title.y = element_text(size = 22),
# #         axis.title.x = element_text(size = 22),
# #         plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
# #         axis.text.x = element_text(size=18),
# #         axis.text.y=element_text(size=16))
# # 
# # p1_leg <-
# #   ggplot() +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Positive droplet"), size = 1) +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Negative droplet"), size = 1) +
# #   theme_bw()+
# #   scale_color_manual(
# #     name='Detection',
# #     breaks=c("Positive droplet", "Negative droplet"),
# #     values = c("forestgreen", "tomato2"),
# #     guide = guide_legend(
# #       override.aes = list(lty = c(1,1),
# #                           size = c(2,2)))) +
# #   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
# #         legend.key.width = unit(1.2,"cm"),
# #         legend.title = element_text(size=18),
# #         legend.key.height = unit(0.8, 'cm'),
# #         legend.text = element_text(size=16))
# # 
# # p2_leg <-
# #   ggplot() +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Cod + herring"), size = 1) +
# #   geom_point(aes(x = NA, y = NA, 
# #                  color = "Cod + saithe"), size = 1) +
# #   theme_bw()+
# #   scale_color_manual(
# #     name="Assay",
# #     breaks=c("Cod + herring", "Cod + saithe"),
# #     values = c("black", "black"),
# #     guide = guide_legend(
# #       override.aes = list(shape = c(17,19),
# #                           size = c(2,2)))) +
# #   theme(plot.margin = margin(1, 4, 1, 0, "cm"),
# #         legend.key.width = unit(1.2,"cm"),
# #         legend.title = element_text(size=18),
# #         legend.key.height = unit(0.8, 'cm'),
# #         legend.text = element_text(size=16))
# # 
# # p1_leg_leg <- cowplot::get_legend(p1_leg+
# #                                     
# #                                     theme(legend.justification = 
# #                                             c(0,0.0)))
# # p2_leg_leg <- cowplot::get_legend(p2_leg+
# #                                     theme(legend.justification = 
# #                                             c(0,1.0)))
# # 
# # legend <- cowplot::plot_grid(p1_leg_leg,
# #                              p2_leg_leg,
# #                              nrow = 2)
# # 
# # Figure_S4 <- cowplot::plot_grid(pp1,
# #                                 legend,ncol = 2, 
# #                                 rel_widths = c(4,1))
# # 
# # Figure_S4
# # 
# # # save the figure if the above 'bSaveFigures' is TRUE
# # if(bSaveFigures==T){
# #   ggsave(plot = Figure_S4, 
# #          # define the output filenmae by pasting together 
# #          # qpcrrunno and qpcrrundate
# #          filename = paste0(outdir02,"/FigS04_Pos_Neg_droplets",
# #                            ".png"),
# #          width=210*1.6,height=297*0.8,
# #          #width=297,height=210,
# #          units="mm",dpi=300)
# # }
# # # 
# # 
# # 
# # #5.5 Figure S5 - LoD and LoQ
# # 
# # ## Filter out 2 outliers
# # DAT <- st_qpcr %>% 
# #   dplyr::filter(!(Well=='B2'&Ct==16.80747223)) %>% 
# #   dplyr::filter(!(Well=='G1'&Ct==44.88344574))
# # 
# # #DAT %>% dplyr::rename("Ct"="Cq")
# # 
# # #View(DAT)
# # pp5 <- DAT %>%
# #   ggplot() +
# #   geom_point(aes(x = log10(int_concentation), 
# #                  y = Ct, color = Species), size = 2) +
# #   facet_wrap(~Species) +
# #   scale_color_manual(values = 
# #                        c("tomato2", "deepskyblue2","orange2")) +
# #   scale_shape_manual(values = 
# #                        c(16, 17, 18)) +
# #   theme_bw()
# # pp5  
# # #source(here('Code','Surpressed','Figure_S5_LoD_LoQ.R'))
# # 
# # #Figure_S5