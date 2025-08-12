#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# R-code  for :
# “Infering opitmal primer concentrations from qPCR MxPro xls files”
# Authors: Steen Wilhelm Knudsen.

# # to use ggplot
# if(!require(ggplot2)){
#   install.packages("ggplot2")
# }
library(ggplot2)

# #library(readxl) # to read the excel file
# if(!require(readxl)){
#   install.packages("readxl")
# }
library(readxl)

# #library(tidyr) # for the seperate function
# if(!require(tidyr)){
#   install.packages("tidyr")
# }
library(tidyr)
# if(!require(readr)){
#   install.packages("readr")
# }
library(readr)
# #install.packages("pals") # to use kovesi.rainbow function
# if(!require(pals)){
#   install.packages("pals")
# }
library(pals)

#library(pals) # to use kovesi.rainbow function

#get the plyr package to run the function 'desc'
# #install.packages("plyr")
# if(!require(plyr)){
#   install.packages("plyr")
# }
library(plyr)

#remove everything in the working environment, without a warning!!
#rm(list=ls())
#Make a function that uses the Convert Integer Vectors to or from UTF-8-encoded Character Vectors
letter2num <- function(x) {utf8ToInt(x) - utf8ToInt("a") + 1L+32}
# set working directory
wd00 <- getwd()
#wd00 <- "/Users/steenknudsen/Documents/Documents/MS_amphibian_eDNA_assays/MS_suppm_amphibia_eDNA"
#wd00 <- "/home/hal9000/Documents/shrfldubuntu18/compare_eDNA_in_ddPCR_and_qPCR"

# setwd (wd00)
# getwd()
#define output directory
wd02.1 <- "output02_plots_from_primer_optimal_conc"
#define directory with input flies
#wd01 <- "inp02_primer_opt_xls_amphibia"
wd01 <- "data/data_qpcr_runs"
# define full path for input directory
#inpfdir01 <- wd00
inpfdir01 <- paste(wd00,"/",wd01,sep="")
# define an outout directory

#wd03 <- paste(wd00,"/",wd02.1,"/",wd02.2,sep="")
wd03 <- paste(wd00,"/",wd02.1,sep="")
#delete previous versions of the output directory
unlink(wd03, recursive=TRUE)
#create new directory
dir.create(wd03)
#make a list of the input csv files
csv.qpcr.fls <- list.files(path=inpfdir01, 
                           pattern="*.csv", full.names=TRUE, recursive=FALSE)
#make a list of the input csv files
txtrep.qpcr.fls <- list.files(path=inpfdir01, 
                           pattern="*Text Report Data", full.names=TRUE, recursive=FALSE)


#make a list of the input xls files
xls.qpcr.fls <- list.files(path=inpfdir01, 
                           pattern="*.xls", full.names=TRUE, recursive=FALSE)
#s
txtrep.qpcr.fls <- txtrep.qpcr.fls[grepl("primer_optim_",txtrep.qpcr.fls)]
#xls.qpcr.fls <- xls.qpcr.fls[grepl("Bombom",xls.qpcr.fls)]

#make an empty list to hold species, qpcr no and plots
ls_sp <- ls()
ls_qpcr <- ls()
ls_plotno <- ls()
#make a starting number to add to when iterating over files
i <- 1
#iterate over files
#for (inf01 in xls.qpcr.fls)
for (inf01 in txtrep.qpcr.fls)
#{print( inf01)}
{
#define one input file
#inf01 <- xls.qpcr.fls[2]
# read in input file
df1<- read.table(inf01, header = T, sep="\t")
# df1<- read_excel(inf01)
# df1 <- as.data.frame(df1)
# 
df1$WellName2 <- df1$Well_1
df1$WellName2 <- df1$Well.Name
#replace in column if there are numbers on species abbreviations
df1$WellName2 <- gsub("(^[A-Za-z]{6})[0-9]+.*_([A-Za-z].*uM.*uM.*uM.*$)","\\1_\\2",df1$WellName2)

#df1$WellName2
#get PCR no and species naem from file
qPCRno <- gsub("^.*_qpcr([0-9]+)_.*$","\\1",gsub("*.*/(.*$)","\\1",inf01))
qPCRno <- gsub("^qpcr([0-9]+)_.*$","\\1",gsub("*.*/(.*$)","\\1",inf01))
specnm <- gsub("^qpcr([0-9]+)_prim_optim_([A-Za-z]{6})_.*$","\\2",gsub("*.*/(.*$)","\\1",inf01))
specnm <- gsub("qpcr[0-9]+_([A-Za-z]{6})_prim.*$","\\1",specnm)

# define assay name and qPCR number
qPCRNospec <- paste("qpcr",qPCRno,"_",specnm,sep="")
#replace spaces in column names
colnames(df1) <- gsub(" ","",colnames(df1))
#new variable to the original dataframe
df1$WellName <- df1$Well_1
df2 <-df1

#separate 
df1<-df1 %>%
  tidyr::separate("Well.Name", into = c("Well.Name","Quantity"), sep=" ")
df1<-df1 %>%
  tidyr::separate(Well.Name, into = c("species","FWN1", "RWN2","PWN3"), sep="_")

#check if the primers have been listed in reverse order
if (unique(substr(df1$FWN1, 1, 1))=="R")
{
##  comment this section out if order of R and F primer is as in df1 - start comment
##  separate differently , if R and F primer were added differently
df2<-df2 %>%
 separate("WellName", into = c("WellName","Quantity"), sep=" ")
df1<-df2 %>%
 separate(WellName, into = c("species","RWN2", "FWN1","PWN3"), sep="_")
## comment this section out if order of R and F primer is as in df1 - end comment
}

#replace regex , replace F with nothing, and uM with nothing
df1$FWN1<-gsub("F","",df1$FWN1)
df1$FWN1<-gsub("uM","",df1$FWN1)
#make variable numeric
df1$FWN1<-as.numeric(df1$FWN1)
#replace regex , replace R with nothing, and uM with nothing
df1$RWN2<-gsub("R","",df1$RWN2)
df1$RWN2<-gsub("uM","",df1$RWN2)
#make variable numeric
df1$RWN2<-as.numeric(df1$RWN2)
#replace regex , replace P with nothing, and uM with nothing
df1$PWN3<-gsub("P","",df1$PWN3)
df1$PWN3<-gsub("uM","",df1$PWN3)
#make variable numeric
df1$PWN3<-as.numeric(df1$PWN3)
#rename specific column names
names(df1)[names(df1) == 'FWN1'] <- 'Forward.concentration'
names(df1)[names(df1) == 'RWN2'] <- 'Reverse.concentration'
names(df1)[names(df1) == 'PWN3'] <- 'Probe.concentration'
#get unique and sort
Forward<-sort(c(unique(df1$`Forward.concentration`)),decreasing = TRUE)
Reverse<-sort(c(unique(df1$`Reverse.concentration`),unique(df1$`Reverse.concentration`)))
df1<-cbind(df1,Row=substr(df1$Well, 1, 1))
df1<-cbind(df1,Column=substr(df1$Well, 2, 3))

df1$Column<-as.numeric(as.character(df1$Column))
#replace in column name
colnames(df1) <- gsub("Ct\\(dRn\\)","Ct..dRn.",colnames(df1))
df1$Ct..dRn. <- df1$Ct
#df1$`Ct (dRn)`<-as.numeric(df1$`Ct (dRn)`)
df1$Ct..dRn.<-as.numeric(df1$Ct..dRn.)


df1$`Ct (dRn)` <- df1$Ct..dRn.
#replace NAs with 1
df1[is.na(df1)] <- 1

df1<-df1[order(df1$Row, df1$Column),]
rownames(df1) <- NULL
###
df1$rownumber<-unname(sapply(as.character(df1$Row), letter2num))
#trasnform in into characters
df1$rowcharacter<-as.character(df1$rownumber)
df1$columncharacter<-as.character(df1$Column)
df1<-df1[order(df1$rownumber, df1$Column),]
#delete rows where `Ct (dRn)` is 1, i.e. remove the wells that did not amplify
df1<-df1[!(df1$`Ct (dRn)`==1),]
#nrow(df1)
#get mean Ct values
df2<-aggregate(df1$"Ct (dRn)",
               by = list(df1$"Forward.concentration",
                         df1$"Reverse.concentration"),
               FUN="mean")
#edit column names
colnames(df2)<-c("Forward.concentration","Reverse.concentration","Ct.dRn")
#get lmit values
max.ct <- max(df2$Ct.dRn)
upppl1 <- ceiling(max.ct)#+0.6
upppl1 <- round(max.ct,digits = 1)+0.2
#define limits
min.ct <- min(df2$Ct.dRn)
lowmi1 <- round(min.ct,digits = 1)-0.2
#lowmi1 <- floor(min.ct)#-0.6
#head(df2,5)
#make a plot
# h<-ggplot(df2, aes(x = reorder(Reverse.concentration,
#                                desc(Reverse.concentration)),
#                    y= reorder(Forward.concentration,
#                               desc(Forward.concentration)))) +

  h<-ggplot(df2, aes(x = reorder(Reverse.concentration,
                                 (Reverse.concentration)),
                     y= reorder(Forward.concentration,
                                (Forward.concentration)))) +
  
  
  geom_tile(aes(fill = Ct.dRn),height=1) +
  geom_text(size=7.0, aes(label = base::round(Ct.dRn, 2),fontface = "bold")) +
  
  # see  https://www.rdocumentation.org/packages/pals/versions/1.6/topics/kovesi
  
  scale_fill_gradientn(colours=kovesi.linear_bgyw_15_100_c68(600), 
                       guide = "colourbar",limits = c(lowmi1,upppl1),
                       values=c(0,1),na.value = "grey50")+
  scale_y_discrete(position = "left")+
  #ggtitle(qPCRNospec)+
labs(x = expression("Reverse.concentration"~paste(mu, "M")), 
       y=expression("Forward.concentration"~paste(mu, "M")))+
theme(axis.text.y = element_text( size = 18,hjust=0.9),
      axis.text.x = element_text( size = 18,hjust=0.5),
      aspect.ratio =0.5,
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text( size = 18),
      legend.title = element_text(size=18),
      legend.text = element_text(size=18))

  #pad with zeros to two characters
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  ipz <-stringr::str_pad(i, 2, pad = "0")
  
  #assign no to plot
  pltno <- paste("plt",ipz,"_",specnm,sep="")
  assign(pltno,h)
  #add to lists
  ls_sp[[i]] <- specnm 
  ls_qpcr[[i]] <- qPCRno
  ls_plotno[[i]] <- pltno
  #add to icreasing number 
  i <- i+1
# set working directory
#setwd (wd00)
pthout <- wd03
#getwd()
#________________________________________________________________________________
# make an eps file with the plot - start
#________________________________________________________________________________
# # Exporting EPS files via postscript()
# postscript(c(paste(pthout,"plot_optimal_primer_conc_",qPCRNospec,".eps", sep = "")),
#            width=(1.6*8.2677),height=(1.5*1.6*2.9232),
#            #family = "Arial", 
#            paper = "special", onefile = FALSE,
#            horizontal = FALSE)
# 
# print(h)
# # end pdf file to save as
# dev.off()
#________________________________________________________________________________
# make an eps file with the plot - end
#________________________________________________________________________________


pdf(c(paste(pthout,"/heatmap_opt_primer_conc_",qPCRNospec,".pdf",  sep = ""))
    ,width=(1.6*8.2677),height=(1.6*2.9232*2))
print(h)
dev.off()
# end iteration over input files
}

#end any previous plots
# dev.off()
##
#combine lists to a data frame
df_plots01 <- as.data.frame(as.matrix(cbind(ls_qpcr, ls_sp, ls_plotno)))
#change column names
colnames(df_plots01) <- c("qpcrNo","spc","pltNo")
# Sort by vector name [spc] then [qpcrNo]
df_plots02 <- df_plots01[
  with(df_plots01, order(spc, qpcrNo)),
  ]
#make a column numeric
df_plots02$qpcrNo <- as.numeric(df_plots02$qpcrNo)
#colnames(df_plots02)
#https://community.rstudio.com/t/how-to-select-top-n-highest-value-and-create-new-column-with-it/38914
#https://stackoverflow.com/questions/24237399/how-to-select-the-rows-with-maximum-values-in-each-group-with-dplyr
#load the dplyr package
library(dplyr)
#use the dplyr to group by 'spc' and then filter inside each group
tibl_plt03 <- df_plots02 %>% 
  dplyr::group_by(spc) %>%
  #select among the highest qPCRno
  #because the highest number equals the most recent qPCR performed
  dplyr::filter(qpcrNo == max(qpcrNo))

plt.empty <- ggplot() + theme_void()
if (!((nrow(tibl_plt03)) %% 2) == 0)
{
  
  plt.empty <- ggplot() + theme_void()
  df_tpl03 <- as.data.frame(tibl_plt03)
  df_tpl03 <- rbind(df_tpl03,NA)
  df_tpl03$pltNo[is.na(df_tpl03$pltNo)] <- plt.empty
  tibl_plt03 <- as_tibble(df_tpl03)
}
tibl_plt03$pltNo[is.na(tibl_plt03$qpcrNo)] <- plt.empty
#make the tibble a data frame
df_plt04 <- as.data.frame(tibl_plt03)

# arrange in columns with two elements per column
# because it fits nicely with two plots in one column on an A4 page
df_plt04 <-  as.data.frame(matrix(unlist(df_plt04$pltNo), nrow=2))
#https://stackoverflow.com/questions/25401111/left-adjust-title-in-ggplot2-or-absolute-position-for-ggtitle
library(ggplot2)
library(grid)
library(gridExtra)
# prepare titles to use for subfigure letters in grid arranged plots
#make a title for one plot to use in the grid arrange
title.grob01 <- textGrob(
  label = "A",
  x = unit(0.2, "lines"), 
  y = unit(0.2, "lines"),
  hjust = 0, vjust = 0,
  gp = gpar(fontsize = 32))
#make a title for the second plot to use in the grid arrange
title.grob02 <- textGrob(
  label = "B",
  x = unit(0.2, "lines"), 
  y = unit(0.2, "lines"),
  hjust = 0, vjust = 0,
  gp = gpar(fontsize = 32))
# get the number of columns from the df with plots
nclpl <- ncol(df_plt04)

#iterate over a sequence
for (i in seq(1:nclpl))
{
  print(i)
#get the two plots per column
plt01 <- get(as.character(df_plt04[,i][1]))
plt02 <- get(as.character(df_plt04[,i][2]))
#arrange plots in grid
p3 <-grid.arrange(arrangeGrob(plt01, top = title.grob01),
                  arrangeGrob(plt02, top = title.grob02),
                  #top = "Global Title", ncol=1) #Use this if you want a global title
                  top = " ", ncol=1)
# https://stackoverflow.com/questions/25401111/left-adjust-title-in-ggplot2-or-absolute-position-for-ggtitle
# https://statisticsglobe.com/arrange-list-of-ggplot2-plots-in-r
# https://stackoverflow.com/questions/10581440/error-in-grid-calll-textbounds-as-graphicsannotxlabel-xx-xy-polygon
# https://stackoverflow.com/questions/53340828/add-title-using-grid-arrange-for-multiple-plots-made-with-gridextragrid-arrang

#define an output file
# pdf(c(paste(pthout,"/heatmap_opt_primer_conc_0",i,".pdf",  sep = ""))
#     ,width=(1.6*8.2677),height=(1.6*2.9232*2))
outfl <- paste(pthout,"/heatmap_opt_primer_conc_0",i,".png",  sep = "")
# png(c(paste(pthout,"/heatmap_opt_primer_conc_0",i,".png",  sep = ""))
#     ,width=(1.6*210),height=(1.2*297),
#     res=150,
#     units = "mm")

ggsave(plot=p3, 
       filename=outfl,width=(1.6*210),height=(1.2*297),
       units = "mm")
# 
# #draw the plot
# grid::grid.draw(p3)
# #close the pdf
# dev.off()
#end iterating over columns in data frame with plots
}


#
#
#