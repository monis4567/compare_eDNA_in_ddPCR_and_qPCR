#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code  for :
# 
#
# “Standard dilution curves from qPCR MxPro txt files”
#
# Authors: Steen Wilhelm Knudsen.

#

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
#
# must be located in the same working directory - as specified in the code below
#
#This code is able to run in:
#
#R studio: Version 0.98.994 – © 2009-2013 RStudio, Inc.
#Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_3) AppleWebKit/537.76.4 (KHTML, like Gecko)
#
#
#____________________________________________________________________________#

#see this
#website
#on how to only install required packages
# #https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(
#   scales, 
#   fields, 
#   gplots,
#   plyr)

library(plyr)
library(scales)
library(gplots)
library(fields)
library(dplyr)



## install the package 'scales', which will allow you to make points on your plot more transparent
#install.packages("scales")
# if(!require(scales)){
#   install.packages("scales")
# }
library(scales)

# get library for filling strip headers on facet wrap ggplots
# https://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set
# if(!require(ggh4x)){
#   install.packages("ggh4x")
# }
library(ggh4x)

#install.packages("fields")
# if(!require(fields)){
#   install.packages("fields")
# }
library(fields)

## install the package 'gplots', to be able to translate colors to hex - function: col2hex
#install.packages("gplots")
# if(!require(gplots)){
#   install.packages("gplots")
# }
library(gplots)

## install the package 'glad', to be able to color using the function 'myPalette'
#install.packages("glad")
#library(glad)

require(graphics)

#get package to read excel files
# if(!require(readxl)){
#   install.packages("readxl")
# }
library(readxl)

#get package to do count number of observations that have the same value at earlier records:
# see this website: https://stackoverflow.com/questions/11957205/how-can-i-derive-a-variable-in-r-showing-the-number-of-observations-that-have-th
#install.packages("plyr")
# if(!require(plyr)){
#   install.packages("plyr")
# }
library(plyr)


##########################################################################################
# begin -  Function to fill NAs with previous value
##########################################################################################
#fill NAs with latest non-NA value
#http://www.cookbook-r.com/Manipulating_data/Filling_in_NAs_with_last_non-NA_value/
#https://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value

fillNAgaps <- function(x, firstBack=FALSE) {
  ## NA's in a vector or factor are replaced with last non-NA values
  ## If firstBack is TRUE, it will fill in leading NA's with the first
  ## non-NA value. If FALSE, it will not change leading NA's.
  
  # If it's a factor, store the level labels and convert to integer
  lvls <- NULL
  if (is.factor(x)) {
    lvls <- levels(x)
    x    <- as.integer(x)
  }
  
  goodIdx <- !is.na(x)
  
  # These are the non-NA values from x only
  # Add a leading NA or take the first good value, depending on firstBack   
  if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
  else             goodVals <- c(NA,            x[goodIdx])
  
  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1
  
  x <- goodVals[fillIdx]
  
  # If it was originally a factor, convert it back
  if (!is.null(lvls)) {
    x <- factor(x, levels=seq_along(lvls), labels=lvls)
  }
  
  x
}
##########################################################################################
# end -  Function to fill NAs with previous value
##########################################################################################


# # set working directory
# wd00 <- "/Users/steenknudsen/Documents/Documents/MS_amphibian_eDNA_assays/MS_suppm_amphibia_eDNA"
# wd00 <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR"
wd00 <- getwd()
# define dir w output files
wd03 <- "data/MONIS6_2021_data"

wd00_wd03 <- paste0(wd00,"/",wd03)

wdin01.1 <- "data/data_qpcr_runs"
#define directory with input flies
wdin01.2 <- "inp01_speci_ampl_plots_amphibia"
# define full path for input directory
inpfdir01 <- paste(wd00,"/",wdin01.1,sep="")
# define an outout directory
wdout02.2 <- "output05_specificity_amplific_plots"
wdout02.2 <- paste(wd00,"/",wdout02.2,sep="")
#delete previous versions of the output directory
unlink(wdout02.2, recursive=TRUE)
#create new directory
dir.create(wdout02.2)
#make a list of the input files
xls.qpcr.fls <- list.files(path=inpfdir01, 
                           pattern="*.xls", full.names=TRUE, recursive=FALSE)
# grep for only the flies that include the 'spectest' in the name 
xls.qpcr.fls <- xls.qpcr.fls[grepl("spectest",xls.qpcr.fls)]
# How to add to a matrix initiated before the loop starts
# see https://stackoverflow.com/questions/13442461/populating-a-data-frame-in-r-in-a-loop
#number of files
ntrf <-  length(xls.qpcr.fls)
#prepare a matrix that matches
mtrx_fls01 <- matrix(ncol=2, nrow=ntrf)
#make a variable with the element you want to search for
id1 <- "xls"
#grep for this variable in the list -  see this example: https://stackoverflow.com/questions/35880242/r-selecting-element-from-list
ls.fl01.xls <- xls.qpcr.fls[grep(paste0(id1), xls.qpcr.fls)]
#get the number of elements in the list
nos.xls.fls <- length(grep(".xls", xls.qpcr.fls))
#make a sequence
nos.in.ls <- seq(1:nos.xls.fls)
#combine to a dataframe
list.xls.inp01 <- as.data.frame(cbind(xls.qpcr.fls, nos.in.ls))
#put one column from this dataframe in to a list, to be used to loop over
#files <- mxpro_ampl.plot.files$filenm1
files <- xls.qpcr.fls
files <- as.vector(files)
#
####################################################################################
#Function for repating rows
#from this website:
#https://www.r-bloggers.com/a-quick-way-to-do-row-repeat-and-col-repeat-rep-row-rep-col/
####################################################################################
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
####################################################################################
#define input path and file for xls file with assays
inf06 <- paste0(wd00_wd03,"/list_of_specific_assays_MONIS6.xlsx")
df_ass01 <- readxl::read_xlsx(inf06)
# narrow down to 'Krebsdyr'
df_ass01 <-  df_ass01[grepl("Krebs",df_ass01$Organismegruppetype),]
# define the genera names  for the new assays 
NwCrucAss <- c("Homarus",
               "Paralithodes",
               "Callinectes",
               "Hemigrapsus")
# collapse this vector
gen_to_incl <- paste(NwCrucAss, collapse = "|")
# and use it to grep in the column with genera
df_ass01 <-  df_ass01[grepl(gen_to_incl,df_ass01$Genus),]
# exclude assays without a primer
df_ass01 <- df_ass01[!is.na(df_ass01$FprimNm),]
# make a Latin species name with an undescore
df_ass01$Latinsk_navn_w.uSc<- gsub(" ","_",df_ass01$Latinsk_navn)
#define species names
spcls <- df_ass01$Latinsk_navn_w.uSc
#use gsub to split between underscores
gnnm <- gsub("(.*)_.*","\\1",spcls) 
spnm <- gsub("(.*)_(.*)","\\2",spcls)
#combine to a dataframe
df_spcabbr01 <- as.data.frame(cbind(gnnm,spnm))
#get only the first 3 letters
df_spcabbr01$abbrgnnm <- substr(df_spcabbr01$gnnm,1,3)
df_spcabbr01$abbrspnm <- substr(df_spcabbr01$spnm,1,3)
#paste together
df_spcabbr01$Abbrnm <- paste(df_spcabbr01$abbrgnnm,df_spcabbr01$abbrspnm,sep="")
df_spcabbr01$Long_Lt_nm <- paste(df_spcabbr01$gnnm,"_",df_spcabbr01$spnm,sep="")
#get only the first letter
df_spcabbr01$flgnnm <- substr(df_spcabbr01$gnnm,1,1)
#make a short spc nm
df_spcabbr01$shspcnm <- paste(df_spcabbr01$flgnnm,". ",df_spcabbr01$spnm,sep="")

#get a single file from a list
files <- xls.qpcr.fls[4]
files <- xls.qpcr.fls

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# make a function that can make the text in the legend in italics
# https://stackoverflow.com/questions/59554096/ggplot2-italics-in-the-legend
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(x == plain, "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#
##____________________________________________________________________________________
# start- loop over filenames in dataframe prepared from list above
##____________________________________________________________________________________
#make  a number to ad to while terating over files
i <- 1
#make empty lists
lqnos <- list()
lspc <- list()
lplotno <- list()
linpf  <- list()
ldf2 <- list()
#iterate over files
for (mxpro_ampl.plot.filename in files){
  inpfilnm <- gsub("*.*/(*.*)","\\1",mxpro_ampl.plot.filename)
  # get the index number for the filename
  idxnF<- which(mxpro_ampl.plot.filename==files)
  #get qpcrno
  qpcrno <- gsub("^*.*qpcr([0-9]+).*.$","\\1",inpfilnm)
  #get species abbreviation in file name
  flnmspab <- gsub("^*.*qpcr([0-9]+)_([A-Za-z]{6}).*.$","\\2",inpfilnm)
  #match the long species name to the abbreviated name
  long_spc_mn <- df_spcabbr01$Long_Lt_nm[match(flnmspab,df_spcabbr01$Abbrnm)]
  # replace underscores with spaces
  long_spc_mn2 <-  gsub("_"," ",long_spc_mn)
  #get a long title for the plot
  plt_titl_w_long_spc_mn <- paste("qpcrno",qpcrno,"targeting",long_spc_mn,sep=" ")
  #make a list of the input files
  plstp.xls.qpcr.fls <- list.files(path=inpfdir01, 
                             pattern="plate_setup.*.xls", full.names=TRUE, recursive=FALSE)
  xlssht <- plstp.xls.qpcr.fls[grepl(qpcrno,plstp.xls.qpcr.fls)]
  
  library(readxl)
  spsht.stup <- readxl::read_excel(xlssht)
  # get the row number where the setup of well starts
  rwn.wtbs <- which("Well"==spsht.stup[,1])
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
  df_dds$WellNo <- df_dds$Well
  # identify NA named columns and remove
  df_dds <- df_dds[!is.na(colnames(df_dds))]
  df_dds<- df_dds[!grepl("^NA",colnames(df_dds))]
  
  # # read in xls files as tibble
  tib01 <- read_excel(mxpro_ampl.plot.filename)
  #change to a dataframe
  df01 <- as.data.frame(tib01)
  #delete row number 1 -  as this only contains 'NA'
  #https://stackoverflow.com/questions/7942519/deleting-every-n-th-row-in-a-dataframe
  df02 <- df01[-1,]
  #the data frame is off, and need the first column value to be the second value
  df02[2,1] <- df02[1,1]
  #replace the first column value with a column heading
  df02[1,1] <- "well"
  #chnage column names
  colnames(df02)<- df02[1,]
  #put the df back in to the same name but without the first row
  df02 <- df02[-1,]
  # use the function defined above
  #Fill NA gaps on a selected column
  df02$well <- fillNAgaps(df02$well)
  #delete row if the column 'Cycles' matches the word 'Cycles'
  df02 <- df02[!(df02$Cycles=="Cycles"),]
  #convert to numeric
  df02$Cycles  <- as.numeric(df02$Cycles)
  # Rename column where names is "Fluorescence (dRn)"
  names(df02)[names(df02) == "Fluorescence (dRn)"] <- "Fluorescence_dRn"
  #convert to numeric
  df02$Fluorescence_dRn <- as.numeric(df02$Fluorescence_dRn)
  #replace Qty and up until comma with space
  df02$well <- gsub(" Qty*.*,"," ",df02$well)
  #exclude rows where 'ROX' appear
  df02 <- df02[!grepl("ROX",df02$well),]
  df02$well <- gsub(" FAM",", FAM ",df02$well)
  df02$well <- gsub(",,",",",df02$well)
  #unique(df02$well)
  well.splt01 <- data.frame(do.call('rbind', strsplit(as.character(df02$well),',',fixed=TRUE)))
  #nrow(well.splt02)
  well.splt02 <- data.frame(do.call('rbind', strsplit(as.character(well.splt01$X3),'_',fixed=TRUE)))
  well.splt03 <- data.frame(do.call('rbind', strsplit(as.character(well.splt01$X2),' ',fixed=TRUE)))
  well.splt03 <- data.frame(do.call('rbind', strsplit(as.character(well.splt01$X3),' ',fixed=TRUE)))
  well.splt04 <- data.frame(do.call('rbind', strsplit(as.character(well.splt03$X2),'_',fixed=TRUE)))
  #append back to dataframe
  #also replace using 'sub' and convert to factor
  df02$repl.no <- as.numeric(as.factor(sub("Repl. ", "",well.splt01$X2)))
  #check if the xls-plot-file is 'wrong' i.e. without any rows, and needs to be skipped
  if(length(df02$repl.no)==0)
  {print("repl.no is 0")} else {
    #}
    #check 
    if (is.null(well.splt01$X4)) {
      df02$probe.col <- as.factor(well.splt01$X3)} else {
        df02$probe.col <- well.splt01$X4}
    
    #append more back to the original dataframe
    if (is.null(well.splt02$X1)) {
      df02$repl.symb <- as.factor(well.splt01$X1)} else {
        df02$repl.symb <- well.splt02$X1}
    if (is.null(well.splt02$X1)) {
      df02$spc.abbr <- as.factor(well.splt01$X2)} else {
        df02$spc.abbr <- well.splt02$X2}
    df02$well.type <- well.splt01$X2
    df02$well.vol <- well.splt02$X4
    df02$well.no <- well.splt01$X1
    df02$spctest_flcol <- well.splt01$X2
    # check the spcs name is in the df of abbreviated names, 
    #if yes, then use the long name, if not then do not change
    df02$spc.abbr <- ifelse(df02$spc.abbr %in% df_spcabbr01$Abbrnm,
                            df_spcabbr01$shspcnm[match(df02$spctest,
                                df_spcabbr01$Abbrnm)],
                            gsub(" ","",as.character(well.splt02$X1)))
    
    df02$repl.symb <-  df02$spc.abbr
    well.splt05 <- data.frame(do.call('rbind', 
                                      strsplit(as.character(df02$spctest_flcol)
                                               ,' ',fixed=TRUE)))
    
    #check if the column has 'Repl'
    df02$spctest <- ifelse(grepl("Repl",df02$spctest),
                           as.character(df02$repl.symb),
                           as.character(well.splt05$X2))
    df02$spctest <- gsub(" ","",df02$spctest)
    # #replace numbers
    df02$spctest <- (gsub("[0-9]+","",df02$spctest))
    # #replace unneeded spaces
    df02$spctest <- (gsub(" ","",df02$spctest))
    # check the spcs name is in the df of abbreviated names, if yes, then use the long name, if not then do not change
    df02$spctest2 <- ifelse(df02$spctest %in% 
                              df_spcabbr01$Abbrnm, 
                            df_spcabbr01$shspcnm[match(df02$spctest,
                                                       df_spcabbr01$Abbrnm)],
                            df02$spctest)
    #check if the string needs to be split
    df02$spctest2 <- ifelse(grepl("_",df02$spctest2),
                            as.character(data.frame(do.call('rbind',
                                strsplit(as.character(df02$spctest2)
                                         ,'_',fixed=TRUE)))[,1]),
                            df02$spctest2)
    #check again
    df02$spctest2 <- ifelse(df02$spctest2 %in% 
                              df_spcabbr01$Abbrnm, 
                            df_spcabbr01$shspcnm[match(df02$spctest2,
                            df_spcabbr01$Abbrnm)],
                            df02$spctest2)
    #check if the species is NA
    if (any(is.na(df02$spctest2))) {df02$spctest2 <- well.splt02$X2}
    
    #unique(df02$spctest2)
    #replace unneeded spaces
    df02$spctest2 <-gsub(" ","",df02$spctest2)
    if (any(grepl("Repl",df02$spctest2))) {df02$spctest2 <- df02$spc.abbr}
    # check if any names are NA, and use a different column in case they are
    if (any(is.na(df02$spc.abbr))) {df02$spc.abbr <- df02$spctest}
    
    
    #View(df02)
    #get a column for a colur
    df02$fl_col <- well.splt05$X5
    #remove rows that has 'ROX'
    df02 <- df02[!grepl("ROX",df02$well),]
    df02$well.conc1 <- df02$well.type
    #replace to w gsub to get only file name
    inpflnm01 <- gsub("*.*/(*.*$)","\\1",mxpro_ampl.plot.filename)
    # make a common column to be able to use dplyr left_join
    df02$WellNo <- df02$well.no
    # use dplyr left_join
    library(dplyr)
    df02 <- df02 %>% dplyr::left_join(df_dds,by="WellNo")
    colnames(df02) <- gsub(" ","_",colnames(df02)) 
    # substitute in long names w content
    fpNm <- df02$ReplNo_Spec_dil_vol_primerprobe
    fpNm <- gsub("(^[A-Za-z]{+}[0-9]{+})_.*","\\1",fpNm)
    fpNm <- gsub("(^[A-Za-z]{+})_.*","\\1",fpNm)
    df02$spcNm3 <- fpNm
    
    # if any spec abbrevations have an NA, then use the spcNm3 instead
    if(any(is.na(df02$spc.abbr))) {df02$spc.abbr <- df02$spcNm3}
    # copy into a new column
    df02$spc.abbr4 <- df02$spc.abbr
    # make it a vector for now to replace in using gsub
    mmn <- df02$spc.abbr4
    #Retain only up to 1st occurunce with gsub
    #https://stackoverflow.com/questions/13608988/grab-from-beginning-to-first-occurrence-of-character-with-gsub
    mmn <- gsub("(^.*?)(_.*?)(_.*?).*", "\\1\\2", mmn)
    # substitute other names as well
    mmn[grepl("E_sinensis",mmn)] <- "E_sinensis"
    mmn[grepl("Erisin",mmn)] <- "E_sinensis" # Eriochir_sinensis
    mmn[grepl("Canpag",mmn)] <- "C_pagurus" # Cancer_pagurus
    mmn[grepl("Nepnor",mmn)] <- "N_norvegicus" # Nephrops_norvegicus
    mmn[grepl("Calsap",mmn)] <- "C_sapidus" # Callinectes_sapidus
    mmn[grepl("Hemsan",mmn)] <- "H_sanguineus" # Hemigrapsus_sanguineus 
    mmn[grepl("Hemtak",mmn)] <- "H_takanoi" # Hemigrapsus_takanoi 
    mmn[grepl("Homame",mmn)] <- "H_americanus" # Homarus_americanus 
    mmn[grepl("Homgam",mmn)] <- "H_gammarus" # Homarus_gammarus 
    mmn[grepl("Carmae",mmn)] <- "C_maenas" # Carcinus maenas 
    mmn[grepl("Hyaara",mmn)] <- "H_araneus"  #Hyas araneus
    mmn[grepl("Rhihar",mmn)] <- "R_harrisii"   # Rhithropanopeus harrisii 
    mmn[grepl("NTC",mmn)] <- "NTC"
    df02$spc.abbr4 <- mmn
    #https://stackoverflow.com/questions/31751022/a-function-to-create-multiple-plots-by-subsets-of-data-frame
    #make a color range for the facet wrap strip headers
    strip.colf <- ggh4x::strip_themed(background_x = elem_list_rect(fill = rainbow(6)))
    
    # make a range of colours for the geom_points in the ggplots
    #http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    # The palette with black:
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                    "#ffa4ff", "#ff209f", "#8c0020", "#06167F",
                    "#B2B23E", "#A56608", "khaki2", "#97FF59",
                    "#065d74")
    # get the unique species names to use in the legend
    # and substitute, so that the names have a punctuation mark
    # added to the abbreviated genus name
    lgndNms <- gsub("_",". ",unique(df02$spc.abbr4))
    # reorder the data frame by the species name
    # using 'dplyr::arrange'  to make sure the plotted lines and species
    # are presented in alphabetical order
    df02 <- df02 %>% dplyr::arrange( spc.abbr4)
    # plot the specificity amplification curves
    library(ggplot2)
    plot01 <- ggplot(df02, aes(x = Cycles,
                               y = Fluorescence_dRn, 
                               group= well.no, 
                               color = spc.abbr4)) +
                               #color = spctest2)) +
      theme_bw() +
      #scale_color_viridis_d(option = 'viridis', direction=-1) +   
      
      #scale_colour_brewer(palette = "Paired") +
      ggplot2::scale_colour_manual(values=cbbPalette,
                                   labels = toexpr(lgndNms, plain = 'Wt')) +
      #https://stackoverflow.com/questions/41631806/change-facet-label-text-and-background-colour/60046113#60046113
      #ggtitle(inpflnm01) +
      ggtitle(long_spc_mn2) +
      # see tip for changing the font on the titel here:
      # https://r-graph-gallery.com/289-control-ggplot2-title.html
      theme(
        plot.title=element_text(family='', face='italic', colour='black', size=12) )+
      #
      geom_point() + 
      # #see this website: https://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set
      # ggh4x::facet_wrap2(~Primer_F~Primer_R~Probe_P , ncol = 2,    
      #            strip = strip.colf) + #'facet_wrap' subsets by column value in dataframe
      # 
      facet_wrap(~Primer_F~Primer_R~Probe_P , ncol = 2) + 
      # # see : https://r-charts.com/ggplot2/facets/
      theme(strip.text = element_text(#face = "bold",
                                      color = "black",
                                      hjust = 0,
                                      size = 8),
            strip.background = element_rect(fill = c("white"),
                                            linetype = "solid",
                                            color = "black",
                                            linewidth = 1)) +
      guides(color=guide_legend(title="species")) +
      labs(y = "Flourescence dRn",
           x= "qPCR cycle") +
      geom_line() #add lines
    
    #pad with zeros to 2 characters for 
    #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
    TwN <- ifelse(nchar(idxnF)<2,stringr::str_pad(idxnF, 2, pad = "0"),idxnF)
    TwN <- as.character(TwN)
    # make a figure number
    FgNo <- paste0("Fig05_",TwN)
    #substitute in the the variable, escape the point with double backslash
    inpflnm02 <- gsub(" ","_",inpflnm01)
    plot.nm2 <- sub("\\.xls", "_",inpflnm02)
    # paste the figure number on to the file name
    plot.nm3 <-  paste(FgNo,"_",plot.nm2,"out02",sep="")
    plot.nm4 <- paste(wdout02.2,"/",plot.nm3,".png",sep="")
    #print the plot in a png
    png(c(plot.nm4)
        ,width=(210),height=(0.6*297),units="mm", res=150)
    print(plot01)
    dev.off()
    pltno <- paste("plot01_no_",i,"_qpcrno",qpcrno,sep="")
    assign(pltno,plot01)
    
    lqnos[[i]] <- qpcrno 
    lspc[[i]] <- flnmspab
    linpf[[i]] <- inpfilnm
    lplotno[[i]] <- pltno
    ldf2[[i]] <- df02
    i <- i+1
    #end if test for repl is 0
  }
  #end iterate over files
}
#_______________________________________________________________________________
# Make a table that list all tissue samples that was used as basis
# for the test of the specificity of the assays
#_______________________________________________________________________________
# define the working directory
wd00 <- getwd()
# define dir w output files
wd03 <- "data/MONIS6_2021_data"
#make complete path to output dir
wd00_wd03 <- paste(wd00,"/",wd03,sep="")
# define input file with list of all assays
inf05 <- "list_of_specific_assays_MONIS6.xlsx"
# paste path and input file together
pthf05 <- paste0(wd00_wd03,"/",inf05)
# read in the file with the list of all assays
tblass <- readxl::read_xlsx(pthf05,col_names = T )
tblass$AbbrvNm
# combine data frames that are stored in the list
# https://stackoverflow.com/questions/2851327/combine-a-list-of-data-frames-into-one-data-frame-by-row
df03  <- bind_rows(ldf2, .id = "column_label")
# only keep unique rows
df04 <- df03  %>% dplyr::distinct(spcNm3,spc.abbr4 )
df04 <- df04 %>% dplyr::arrange(spc.abbr4)
# exclude the NTC
df04 <- df04[!grepl("NTC",df04$spcNm3),]
# make a vector with the abbreviated species names
AbbspcNms <- unique(df04$spc.abbr4)
# make vector with the complete species names
FlspcNms <- c(
  "Carcinus maenus",
  "Cancer pagurus",
  "Callinectes sapidus",
  "Eriocheir sinensis",
  "Homarus americanus",
  "Hyas araneus",
  "Hyas coarctatus",
  "Homarus gammarus",
  "Hemigrapsus sanguineus",
  "Hemigrapsus takanoi",
  "Lithodes maja",
  "Neolithodes asperrimus",
  "Nephrops norvegicus",
  "Pagurus bernhardus",
  "Paralithodes camtschaticus",
  "Rhithropanopeus harrisii"
)
# make a data frame with the abbreviated and full species names
cbnms <- as.data.frame(cbind(AbbspcNms,FlspcNms))
# change the column names in the data frame with tissue samples
colnames(df04) <- c("spcNm3","AbbspcNms")
#  combine the data frames, using left_join
df04 <- df04 %>% dplyr::left_join(cbnms,by="AbbspcNms")


#
#