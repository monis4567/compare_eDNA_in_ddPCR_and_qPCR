#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code  for :
# “Probe optimization curves from qPCR MxPro txt files”


# Change the working directory to a path on your own computer , and run
# the individual parts below to reproduce the diagrams presented in the paper
# All input data required needs to be available as csv-files in the same directory 
# as this R-code use for working directory.
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

# #install.packages("fields")
# if(!require(fields)){
#   install.packages("fields")
# }
library(fields)

library(ggplot2)
# ## install the package 'gplots', to be able to translate colors to hex - function: col2hex
# #install.packages("gplots")
# if(!require(gplots)){
#   install.packages("gplots")
# }
library(gplots)

## install the package 'glad', to be able to color using the function 'myPalette'
#install.packages("glad")
#library(glad)

require(graphics)

# #get package to read excel files
# if(!require(readxl)){
#   install.packages("readxl")
# }
library(readxl)

#get package to do count number of observations that have the same value at earlier records:
# see this website: https://stackoverflow.com/questions/11957205/how-can-i-derive-a-variable-in-r-showing-the-number-of-observations-that-have-th
# #install.packages("plyr")
# if(!require(plyr)){
#   install.packages("plyr")
# }
library(plyr)

library(ggplot2)
library(grid)
library(gridExtra)
# install.packages("ggthemes")
# install.packages("extrafont")
# install.packages("remotes")
library(ggthemes)
library(extrafont)
library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
#extrafont::font_import()
#After that load the fonts:
library(extrafont)
#extrafont::font_import()
# Load fonts
loadfonts(quiet = T)
#fonts()
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


# set working directory
#wd00 <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR"
wd00 <- getwd()
#setwd (wd00)
#getwd()
wd01 <- "data/data_qpcr_runs"
# define full path for input directory
inpfdir01 <- paste(wd00,wd01,sep="/")
# define an outout directory
wd02 <- "/output03_plots_from_probe_optimal_conc"
wd03 <- paste(wd00,wd02,sep="")
#delete previous versions of the output directory
unlink(wd03, recursive=TRUE)
#create new directory
dir.create(wd03)
#make a list of the input files
xls.qpcr.fls <- list.files(path=inpfdir01, 
  pattern="*.xls", full.names=TRUE, recursive=FALSE)

xls.qpcr.fls <- xls.qpcr.fls[grepl("probe",xls.qpcr.fls)]
# How to add to a matrix initiated before the loop starts
# see https://stackoverflow.com/questions/13442461/populating-a-data-frame-in-r-in-a-loop
#number of files
ntrf <-  length(xls.qpcr.fls)

if (!gtools::even(ntrf))
  {
    ntrf <- ntrf+1
  }
#prepare a matrix that matches
mtrx_fls01 <- matrix(ncol=2, nrow=ntrf)
#make a variable with the element you want to search for
id1 <- "xls"
#grep for this variable in the list -  see this example: https://stackoverflow.com/questions/35880242/r-selecting-element-from-list
ls.fl01.xls <- xls.qpcr.fls[grep(paste0(id1), xls.qpcr.fls)]
#get the number of elements in the list
nos.xls.fls <- length(grep(".xls", xls.qpcr.fls))
#nos.xls.fls <- 3
#make a sequence
nos.in.ls <- seq(1:nos.xls.fls)
#nos.in.ls <- seq(1:nos.xls.fls)
#combine to a dataframe
list.xls.inp01 <- as.data.frame(cbind(xls.qpcr.fls, nos.in.ls))
# #Split up a dataframe by number of rows 
# #https://stackoverflow.com/questions/7060272/split-up-a-dataframe-by-number-of-rows
# #it will be put in a list of dataframes
# spl_df01 <- split(list.xls.inp01,rep(1:nos.xls.fls,each=3))
# 
# #take one of these elements from this list and convert to an individual dataframe
# mxpro_ampl.plot.files <- as.data.frame(spl_df01[1])
# #add a column with numbers for plot - here it is from 1 to 4
# mxpro_ampl.plot.files$no.f.plot <-seq(from=1,to=nos.xls.fls,by=1)
# 
# #change the column names
# colnames(mxpro_ampl.plot.files)[1] <- c("filenm1")
# colnames(mxpro_ampl.plot.files)[2] <- c("nos.ls")
#put one column from this dataframe in to a list, to be used to loop over
files <- list.xls.inp01
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


#copy the list
files <- xls.qpcr.fls
files <- files[grepl("probeoptim",files)]
##____________________________________________________________________________________
# start- loop over filenames in dataframe prepared from list above
##____________________________________________________________________________________
i <- 1
lqnos <- list()
lspc <- list()
lplotno <- list()
linpf <- list()
#
files.w.prob.opt <- files[1]
files.w.prob.opt <- files
for (mxpro_ampl.plot.filename in files.w.prob.opt){
  inpfilnm <- gsub("*.*/(*.*)","\\1",mxpro_ampl.plot.filename)
  # print(mxpro_ampl.plot.filename)
  # }
  #get qpcrno
  qpcrno <- gsub("^*.*qpcr([0-9]+).*.$","\\1",inpfilnm)
  #get species abbreviation in file name
  flnmspab <- gsub("^*.*qpcr([0-9]+)_([A-Za-z]{6}).*.$","\\2",inpfilnm)
  #match the long species name to the abbreviated name
  if (exists("df_spcabbr01")){
  long_spc_mn <- df_spcabbr01$Long_Lt_nm[match(flnmspab,df_spcabbr01$Abbrnm)]} else {long_spc_mn <- flnmspab}
  #get a long title for the plot
  plt_titl_w_long_spc_mn <- paste("qpcrno",qpcrno,"targeting",long_spc_mn,sep=" ")
  #match the current filename in the loop with the filename in the dataframe that also holds
#the number for the file - notice that the match-value needs to be converted to a 
#character
# #____
#   no.f.f2 <- mxpro_ampl.plot.files$nos.ls[match(mxpro_ampl.plot.filename, 
#       as.character(mxpro_ampl.plot.files$filenm1))]
# #match between current filename and return the sequence number
# no.fpl1 <- mxpro_ampl.plot.files$no.f.plot[match(mxpro_ampl.plot.filename, 
#       as.character(mxpro_ampl.plot.files$filenm1))]
# #____
# read in xls files as tibble
tib01 <- readxl::read_excel(mxpro_ampl.plot.filename)
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
#plot the curves
#plot(df02$Cycles,df02$Fluorescence_dRn, type="o", col="blue" ) 
#split by delimiter
#https://stackoverflow.com/questions/7069076/split-column-at-delimiter-in-data-frame
#replace Qty and up until comma with space
df02$well <- gsub(" Qty*.*,"," ",df02$well)
well.splt01 <- data.frame(do.call('rbind', strsplit(as.character(df02$well),',',fixed=TRUE)))
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

spc_prob <- data.frame(do.call('rbind',strsplit(as.character(well.splt01$X3)," ",fixed=T)))
spc_prob$X2 <- gsub(" ","",spc_prob$X2)
spc_prob2 <- data.frame(do.call('rbind',strsplit(as.character(spc_prob$X2),"_",fixed=T)))

df02$spc.abbr <- spc_prob2$X1
df02$probe.conc <- as.numeric(gsub("P","",spc_prob2$X2))
df02$probe.conc <- as.character(df02$probe.conc)
data.frame(do.call('rbind',strsplit(as.character(well.splt01$spc_probconc),"_",fixed=T)))
df02$probe.col <- well.splt01$X3
#append more back to the original dataframe
df02$repl.symb <- well.splt02$X1
#df02$spc.abbr <- well.splt02$X2
df02$well.type <- well.splt02$X3
df02$well.vol <- well.splt02$X4
df02$well.no <- well.splt01$X1
df02$spctest_flcol <- well.splt01$X2
well.splt05 <- data.frame(do.call('rbind', strsplit(as.character(df02$spctest_flcol),' ',fixed=TRUE)))
#check if the column has 'Repl'
df02$spctest <- ifelse(grepl("Repl",df02$spctest),as.character(df02$repl.symb),as.character(well.splt05$X2))
# #replace numbers
 df02$spctest <- (gsub("[0-9]+","",df02$spctest))
# #replace unneeded spaces
 df02$spctest <- (gsub(" ","",df02$spctest))
# check the spcs name is in the df of abbreviated names, if yes, then use the long name, if not then do not change
 if (exists("df_spcabbr01")){
    df02$spctest2 <- ifelse(df02$spctest %in% df_spcabbr01$Abbrnm, df_spcabbr01$shspcnm[match(df02$spctest,df_spcabbr01$Abbrnm)],df02$spctest)
    } else {df02$spctest2 <- df02$spctest}
 #check if the string needs to be split
 df02$spctest2 <- ifelse(grepl("_",df02$spctest2),
   as.character(data.frame(do.call('rbind', strsplit(as.character(df02$spctest2),'_',fixed=TRUE)))[,1]),
   df02$spctest2)
 #check again
 if (exists("df_spcabbr01")){
 df02$spctest2 <- ifelse(df02$spctest2 %in% df_spcabbr01$Abbrnm, df_spcabbr01$shspcnm[match(df02$spctest2,df_spcabbr01$Abbrnm)],df02$spctest2)
 }
 #get a column for a colur
df02$fl_col <- well.splt05$X5
#remove rows that has 'ROX'
df02 <- df02[!grepl("ROX",df02$well),]
#append back to df as numeric and replace with 'gsub' function
#df02$well.vol.val <- as.numeric(sub("uL", "",df02$well.vol))
df02$well.conc1 <- df02$well.type
#replace to w gsub to get only file name
inpflnm01 <- gsub("*.*/(*.*$)","\\1",mxpro_ampl.plot.filename)
#inpflnm01
#df_spcabbr01$Abbrnm
#https://stackoverflow.com/questions/31751022/a-function-to-create-multiple-plots-by-subsets-of-data-frame
library(ggplot2)
plot01 <- ggplot(df02, aes(x = Cycles,
                y = Fluorescence_dRn, 
                group= well.no, 
                color = probe.conc)) +
  geom_point() + 
  theme_bw() +
  scale_color_viridis_d(option = 'viridis', direction=-1) +
  #facet_wrap(~spctest, nrow = 3) + #'facet_wrap' subsets by column value in dataframe
  geom_line() + #add lines
  labs(color="probe conc.\n(µm)") +
  #https://stackoverflow.com/questions/41631806/change-facet-label-text-and-background-colour/60046113#60046113
  #ggtitle(inpflnm01)
  ggtitle(plt_titl_w_long_spc_mn)
#substitute in the the variable, escape the point with double backslash
#plot01
inpflnm02 <- gsub(" ","_",inpflnm01)
plot.nm2 <- sub("\\.xls", "_",inpflnm02)
plot.nm3 <-  paste(plot.nm2,"out02",sep="")
plot.nm4 <- paste(wd03,"/",plot.nm3,".png",sep="")
#print the plot in a png
# png(c(plot.nm4)
#     ,width=(1*210),height=(0.35*297),
#     units = "mm",res=150)
# #print(plot01)
# #draw the plot
# grid::grid.draw(plot01)
# dev.off()

ggsave(plot=plot01, 
       filename=plot.nm4, width=(1*210),height=(0.35*297),
       units = "mm")

#pad with zeros to two characters
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
ipz <-stringr::str_pad(i, 2, pad = "0")

pltno <- paste("plt",ipz,"_qpcrno",qpcrno,sep="")
assign(pltno,plot01)

lqnos[[i]] <- qpcrno 
linpf[[i]] <- inpfilnm
lspc[[i]] <- flnmspab
lplotno[[i]] <- pltno
i <- i+1
#end if test for repl is 0
}
#end iterate over files
}

#

#combine lists to a data frame
df_plots01 <- as.data.frame(as.matrix(cbind(lqnos, lspc,linpf,lplotno)))
#change column names
colnames(df_plots01) <- c("qpcrNo","spc","inpfnm","pltno")
#make lists characters and numeric
df_plots01$qpcrNo <- as.numeric(df_plots01$qpcrNo)
df_plots01$spc <- as.character(df_plots01$spc)
df_plots01$inpfnm <- as.character(df_plots01$inpfnm)
df_plots01$pltno <- as.character(df_plots01$pltno)
#https://community.rstudio.com/t/how-to-select-top-n-highest-value-and-create-new-column-with-it/38914
#https://stackoverflow.com/questions/24237399/how-to-select-the-rows-with-maximum-values-in-each-group-with-dplyr
#load the dplyr package
library(dplyr)
#use the dplyr to group by 'spc' and then filter inside each group
tibl_plt03 <- df_plots01 %>% 
  dplyr::group_by(spc) %>%
  #select among the highest qPCRno
  #because the highest number equals the most recent qPCR performed
  dplyr::filter(qpcrNo == max(qpcrNo))
#make the tibble a data frame
df_plt04 <- as.data.frame(tibl_plt03)
# Sort by vector name [query_def] then [Hsp_bit_score]
df_plt05 <- df_plt04[
  with(df_plt04, order(spc)),
  ]
# check if there are an uneven number of plots
nrw.dpl05 <- nrow(df_plt05)
if (gtools::odd(nrw.dpl05))
{
  clnm.df_plt05 <- colnames(df_plt05)
  cln.dp05 <- ncol(df_plt05) 
  ev <- rep(NA,cln.dp05)
  df_plt05 <- rbind(df_plt05,ev)
  epl <- ggplot()
  df_plt05$pltno[is.null(df_plt05$pltno)] <- epl
}
# because it fits nicely with two plots in one column on an A4 page
df_plt05 <-  as.data.frame(matrix(unlist(df_plt05$pltno), nrow=2))

ndf05 <- which(is.na(df_plt05)) 

#df_plt05[ndf05] 

ggplot() + theme_void()

#https://stackoverflow.com/questions/25401111/left-adjust-title-in-ggplot2-or-absolute-position-for-ggtitle
library(ggplot2)
library(grid)
library(gridExtra)
# prepare titles to use for subfigure letters in grid arranged plots
#make a title for one plot to use in the grid arrange
title.grob01 <- textGrob(
  label = "A",
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0,
  gp = gpar(fontsize = 32))
#make a title for the second plot to use in the grid arrange
title.grob02 <- textGrob(
  label = "B",
  x = unit(0, "lines"), 
  y = unit(0, "lines"),
  hjust = 0, vjust = 0,
  gp = gpar(fontsize = 32))

# # set working directory
# setwd (wd00)
# getwd()
pthout <- paste(wd00,wd02,"/",sep="")
# get the number of columns from the df with plots
nclpl <- ncol(df_plt05)
#iterate over a sequence
for (i in seq(1:nclpl))
{
  print(i)
  
  #get the two plots per column
  plt01 <- get(as.character(df_plt05[,i][1]))
  plt02 <- get(as.character(df_plt05[,i][2]))
  
  #arrange plots in grid
  p3 <-grid.arrange(arrangeGrob(plt01, top = title.grob01),
                    arrangeGrob(plt02, top = title.grob02),
                    #top = "Global Title", ncol=1) #Use this if you want a global title
                    top = " ", ncol=1)
  # https://stackoverflow.com/questions/25401111/left-adjust-title-in-ggplot2-or-absolute-position-for-ggtitle
  # https://statisticsglobe.com/arrange-list-of-ggplot2-plots-in-r
  # https://stackoverflow.com/questions/10581440/error-in-grid-calll-textbounds-as-graphicsannotxlabel-xx-xy-polygon
  # https://stackoverflow.com/questions/53340828/add-title-using-grid-arrange-for-multiple-plots-made-with-gridextragrid-arrang
  
  outfl <- paste(pthout,"plt_opt_probe_conc_0",i,".png",  sep = "")
  
  
  
  ggsave(plot=p3, 
         filename=outfl,width=(1*210),height=(0.4*297*2),
         units = "mm")
  # 
  # #define an output file
  # png(c(paste(pthout,"plt_opt_probe_conc_0",i,".png",  sep = ""))
  #     ,width=(1*210),height=(0.4*297*2),
  #     units="mm",
  #     res=150)
  # #draw the plot
  # grid::grid.draw(p3)
  # #close the pdf
  # dev.off()
  #end iterating over columns in data frame with plots
}

#______________________________________________________________________________
# #______________________________________________________________________________
# #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# # Multiple plot function
# #
# # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# # - cols:   Number of columns in layout
# # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# #
# # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# # then plot 1 will go in the upper left, 2 will go in the upper right, and
# # 3 will go all the way across the bottom.
