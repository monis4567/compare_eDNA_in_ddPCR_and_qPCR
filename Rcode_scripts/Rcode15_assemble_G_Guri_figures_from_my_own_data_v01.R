#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code for the project:
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

# define the working directory
wd00 <- getwd()
#Create a directory to put the code obtained from github into
wd00_wdc <- paste0(wd00,"/code")
#Delete any previous versions of the output directory
unlink(wd00_wdc, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd00_wdc)
#set the output directory
outdir01 <- wd00_wdc

#Create a directory to obtain the data frames from
wd00_wd02 <- paste0(wd00,
                    paste0("/output13_",#,no_sbs,
                           "_my_df_modified_to_match",
                           "_G_Guri_setup"
                    ))
#set the output directory
outdir02 <- wd00_wd02
#define data directory
wd00.1 =paste0(wd00,"/data/MONIS6_2021_data")
# read in an xlsx file with all primer assays listed 
# for each species and the abbreviations
inf04 <- "list_of_specific_assays_MONIS6.xlsx"
wd00.1_inf04 <- paste0(wd00.1,"/",inf04 )
#library(xlsx)
df_dtc_asss <- openxlsx::read.xlsx(wd00.1_inf04,1)
# use the data frame with the detection assays to get the full
# length Latin Genus and Species names
df_dtc_asss$GnNm.l1 <- substr(df_dtc_asss$Genus, 1, 1)
df_dtc_asss$shrtNm <- paste0(df_dtc_asss$GnNm.l1,". ",df_dtc_asss$Species)
#read in the tsv tables
df_ddpcr03 <- read.table(paste0(wd00_wd02,"/df_ddpcr03.tsv"),
                         sep = "\t",
                         header = T) 
df_qpcr03 <- read.table(paste0(wd00_wd02,"/df_qpcr03.tsv"),
                        sep = "\t",
                        header = T)
#in the 'df_qpcr03$Species' and the 'df_ddpcr03$Species' replace "Cragig"
# with 'Maggig', 
df_qpcr03$Species <- gsub("Cragig", "Maggig", df_qpcr03$Species)
df_ddpcr03$Species <- gsub("Cragig", "Maggig", df_ddpcr03$Species)
# Use plyr to select the columns 'shrtNm', 'Phyla', 'Class', 'AbbrvNm', 'Species'
# 'Genus', 'Organismegruppetype', 'Latinsk_navn', 'Navn_dansk' from the 'df_dtc_asss'
# data frame
df_dta <- dplyr::select(df_dtc_asss, shrtNm, Phyla, Class, AbbrvNm, Species,
                        Genus, Organismegruppetype, Latinsk_navn, Navn_dansk)
# limit the 'df_dta' data frame to remove rows where there is duplicate
# occurences of 'AbbrvNm'
df_dta <- df_dta[!duplicated(df_dta$AbbrvNm),]
# then use 'left_join' to merge the 'df_dta' data frame with the 'df_ddpcr03' using
# the 'AbbrvNm' in the 'df_dta' data frame and the 'Species' in the 'df_ddpcr03' data frame
# as the common column
df_ddpcr03 <- dplyr::left_join(df_ddpcr03, df_dta, by = c("Species" = "AbbrvNm"))
# repeat for the 'df_qpcr03' data frame
df_qpcr03 <- dplyr::left_join(df_qpcr03, df_dta, by = c("Species" = "AbbrvNm"))


# get the unique species names from the 'df_ddpcr03' data frame
uqp_spc<- unique(df_qpcr03$Species)
udp_spc<- unique(df_ddpcr03$Species)
# combine the two unique species names in one vector
spcabbr.unq <- unique(c(uqp_spc,udp_spc))
# order them alphabetically
spcabbr.unq <- sort(spcabbr.unq)

# make a new output directory
wd15 <- paste0(wd00,"/output15_assembled_figures_from_G_Guri_test_on_own_data")
#Delete any previous versions of the output directory
unlink(wd15, recursive=TRUE)
#Create a directory to put resulting output files in
dir.create(wd15)
#set the output directory
outdir15 <- wd15

# get a list of the files in the 'wd00_wd02' directory that holds the string 'tibl_podd_'
lst_tbp_fls <- list.files(path = wd00_wd02, pattern = "tibl_podd_", full.names = TRUE)
# get a list of the files in the 'wd00_wd02' directory that holds the string "diff_ddPCR_qPCR_quant_prec_"
lst_dqd_fls <- list.files(path = wd00_wd02, pattern = "diff_ddPCR_qPCR_quant_prec_", full.names = TRUE)
# get a list of the files in the 'wd00_wd02' directory that holds the string "diff_Two_step_qPCR_mod_quant_"
lst_dts_fls <- list.files(path = wd00_wd02, pattern = "diff_Two_step_qPCR_mod_quant_", full.names = TRUE)

# merge all the elements in the list 'lst_tbp_fls' to one single data frame
df_tbp <- 
  lst_tbp_fls %>%
  lapply(read.table, sep = ",", header = TRUE) %>%
  bind_rows()
#in the 'df_tbp$Species' replace "Cragig" with 'Maggig',
df_tbp$Species <- gsub("Cragig", "Maggig", df_tbp$Species)
# then use 'left_join' to merge the 'df_dta' data frame with the 'df_tbp' using
# the 'AbbrvNm' in the 'df_dta' data frame and the 'Species' in the 'df_tbp' data frame
# as the common column
df_tbp <- dplyr::left_join(df_tbp, df_dta, by = c("Species" = "AbbrvNm"))
# merge all the elements in the list 'lst_dqd_fls' to one single data frame
df_dqd <- 
  lst_dqd_fls %>%
  lapply(read.table, sep = ",", header = TRUE) %>%
  bind_rows()
#in the 'df_dqd$Species' replace "Cragig" with 'Maggig',
df_dqd$Species <- gsub("Cragig", "Maggig", df_dqd$Species)
# then use 'left_join' to merge the 'df_dta' data frame with the 'df_dqd' using
# the 'AbbrvNm' in the 'df_dta' data frame and the 'Species' in the 'df_dqd' data frame
# as the common column
df_dqd <- dplyr::left_join(df_dqd, df_dta, by = c("Species" = "AbbrvNm"))

# merge all the elements in the list 'lst_dts_fls' to one single data frame
df_dts <- 
  lst_dts_fls %>%
  lapply(read.table, sep = ",", header = TRUE) %>%
  bind_rows()
#in the 'df_tbp$Species' replace "Cragig" with 'Maggig',
df_tbp$Species <- gsub("Cragig", "Maggig", df_tbp$Species)
# then use 'left_join' to merge the 'df_dta' data frame with the 'df_dts' using
# the 'AbbrvNm' in the 'df_dta' data frame and the 'Species' in the 'df_dts' data frame
# as the common column
df_dts <- dplyr::left_join(df_dts, df_dta, by = c("Species" = "AbbrvNm"))
# In the ggplots , later on, a color is required for each species
# start by making a data frame
# get one range of colors
# get another range of colors
clr02 <- palette.colors(palette = "Polychrome")
clr01 <- palette.colors(palette = "Okabe-Ito")
#clr02 <- palette.colors(palette = "Polychrome")
clr02 <- palette.colors(palette = "Alphabet")

# combine the two ranges of colors
clr03 <- c(clr01,clr02)
# Define the colorblind-friendly palette
clr03 <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#999999", "#AA0DFE", "#3283FE", "#85660D",
  "#782AB6", "#565656", "#1C8356", "#16FF32",
  "#F7E1A0"
)
# use the combined range of colors, to match up with the
# number of species in the data 'tibl_podd' data frame
clr04 <- clr03[1:length(unique(spcabbr.unq))]
col.f_spc <- clr04
spcNm <- spcabbr.unq[order(spcabbr.unq)]

# combine in to data frame
df_clr <- as.data.frame(cbind(spcNm,col.f_spc))
# match the abbreviated name 'spcNm' in the 'df_dtc_asss' data frame
# get the shNm appended on to the 'df_clr' data frame
df_clr$shrtNm <- df_dtc_asss$shrtNm[match(df_clr$spcNm,df_dtc_asss$AbbrvNm)]

# Create a named vector of colors from df_clr
#clr04 <- setNames(df_clr$col.f_spc, df_clr$spcNm)

colnames(df_tbp)
# make a ggplot object
pp1 <-
  df_tbp %>%
  ggplot()+
  geom_line(aes(x=10^C,y=prob_dd,color=shrtNm),lty=2,linewidth=0.9)+
  geom_line(aes(x=10^C,y=prob_q,color=shrtNm),lty=1,linewidth=0.9)+
  # geom_point(data=st_qpcr %>% group_by(Species,Sample_name) %>%
  #              summarise(int_concentation=mean(int_concentation),
  #                        pres=mean(pres)),
  #            aes(x=int_concentation,y=pres,color=Species),shape=15,size=3)+
  # geom_point(data=st_ddpcr %>% group_by(Species,Sample_name) %>%
  #              summarise(int_concentation=mean(int_concentation),
  #                        pres=mean(pres)),
  #            aes(x=int_concentation,y=pres,color=Species),shape=19,size=3)+
  # scale_x_log10(labels=scientific_10,lim=c(1e-4,1e6))+
  scale_x_log10(labels=NULL,lim=c(1e-2,1e6),
                breaks=10^c(-2,-1,0,1,2,4,6))+
  ylim(0,1)+
  # scale_color_manual(values=
  #                      
  #                      #c("tomato2", "deepskyblue2","orange2"))+
  #                      # in the original code only 3 species were included
  #                      # this only required 3 colors. This setup has 17 species
  #                      # and need 17 colors
  #                      c(clr04))+
  scale_color_manual(values = setNames(df_clr$col.f_spc, df_clr$shrtNm)) +
  facet_wrap(~shrtNm, nrow=3)+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
        axis.title.y = element_text(size = 19),
        strip.text = element_text(face = "italic"),
        #plot.margin = margin(0.1, 0.1, -0.5, 0.5, "cm"),
        plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
        axis.title.x = element_text(size = 19),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15))+
  ylab("Detection probability")+
  #xlab("C - Nominal concentration (copies/μL)")
  xlab("")
pp1

bSaveFigures <- T
# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = pp1, 
         # define the output filenmae by pasting together 
         filename = paste0(outdir15,
                           "/Fig15_Nominal_DNA_concentration_",
                           ".png"),
         width=210*1.6,height=297*0.6,
         #width=297,height=210,
         units="mm",dpi=300)
}


theta_line <- df_tbp
pp2 <-
  theta_line %>%
  ggplot()+
  geom_line(aes(x=10^C,y=delta,color=shrtNm),
            lty=6,linewidth=0.9)+
  scale_x_log10(labels=scientific_10,
                lim=c(1e-2,1e6),
                breaks=10^c(-2,-1,0,1,2,4,6))+
  ylim(-0.1,0.8)+
  #scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
  #scale_color_manual(values=c(clr04))+
  scale_color_manual(values = setNames(df_clr$col.f_spc, df_clr$shrtNm)) +
  facet_wrap(~Species,nrow=3)+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 19),
        plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
        axis.title.x = element_text(size = 19),
        axis.text.x = element_text(size=15),
        axis.text.y=element_text(size=15))+
  ylab("ddPCR - qPCR \ndetection probability")+
  xlab("Nominal DNA concentration (copies/uL)")

pp2
# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = pp2, 
         # define the output filenmae by pasting together 
         filename = paste0(outdir15,
                           "/Fig16_Nominal_DNA_concentration_",
                           ".png"),
         width=210*1.6,height=297*0.6,
         #width=297,height=210,
         units="mm",dpi=300)
}
# get the number of unique species in the 'df_tbp' data frame
Nspc <- length(unique(df_tbp$Species))
str(df_clr)
# Create a data frame with species and their corresponding colors
legend_data <- data.frame(
  Species = unique(df_tbp$shrtNm),
  Color = df_clr$col.f_spc[match(unique(df_tbp$shrtNm), df_clr$shrtNm)]
)

# Make a ggplot object
p1_leg <- ggplot(legend_data, aes(color = Color)) +
  geom_point(aes(x = NA, y = NA), size = 1) +
  geom_line(aes(x = NA, y = NA), size = 1) +
  theme_bw() +
  scale_color_manual(
    name = 'Assay',
    breaks = legend_data$Species,
    values = setNames(legend_data$Color, legend_data$Species),
    guide = guide_legend(
      override.aes = list(
        lty = rep(1, nrow(legend_data)),
        size = rep(3, nrow(legend_data))
      )
    )
  ) +
  theme(
    plot.margin = margin(1, 4, 1, 0, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(size = 18),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = 16, face = "italic")
  )

p1_leg_leg <- cowplot::get_legend(p1_leg + theme(legend.justification = c(0, -1.3)))



p1_leg_leg
p3_leg <- ggplot() +
  geom_line(aes(x = NA, y = NA, color = "ddPCR sensitivity"), linewidth = 1) +
  geom_line(aes(x = NA, y = NA, color = "qPCR sensitivity"), linewidth = 1) +
  theme_bw()+
  scale_color_manual(
    name="Modelled sensitivity",
    breaks=c("ddPCR sensitivity","qPCR sensitivity"),
    values = c("black","black"),
    guide = guide_legend(
      override.aes = list(lty = c(2,1),
                          size = c(1,2)))) +
  theme(plot.margin = margin(1, 4, 1, 0, "cm"),
        legend.key.width = unit(1.2,"cm"),
        legend.title = element_text(size=18),
        legend.key.height = unit(0.8, 'cm'),
        legend.text = element_text(size=16))
p3_leg_leg <- cowplot::get_legend(p3_leg+
                                    theme(legend.justification = 
                                            c(0,-0.1)))

p2_leg <- ggplot() +
  geom_line(aes(x = NA, y = NA, 
                color = "Sensitivity difference\n(ddPCR - qPCR)"), 
            size = 1) +
  theme_bw()+
  scale_color_manual(
    name="",
    breaks=c("Sensitivity difference\n(ddPCR - qPCR)"),
    values = c("black"),
    guide = guide_legend(
      override.aes = list(lty = c(6),
                          size = c(1)))) +
  theme(plot.margin = margin(1, 4, 1, 0, "cm"),
        legend.key.width = unit(1.2,"cm"),
        legend.title = element_text(size=1),
        legend.key.height = unit(0.8, 'cm'),
        legend.text = element_text(size=16))
p2_leg_leg <- cowplot::get_legend(p2_leg+
                                    # the 'theme(legend.justification' appears
                                    # to be able to shift the legends on side more around
                                    # it is here changed from 0.95 to 10.95
                                    theme(legend.justification = 
                                            c(0,10.95)))

p4_leg <-
  ggplot() +
  geom_point(aes(x = NA, y = NA, color = "ddPCR standard samples"), size = 1) +
  geom_point(aes(x = NA, y = NA, color = "qPCR standard samples"), size = 1) +
  theme_bw()+
  scale_color_manual(
    name="Measured sensitivity",
    breaks=c("ddPCR standard samples","qPCR standard samples"),
    values = c("black","black"),
    guide = guide_legend(
      override.aes = list(shape = c(19,15),
                          size = c(3)))) +
  theme(plot.margin = margin(1, 4, 1, 0, "cm"),
        legend.key.width = unit(1.2,"cm"),
        legend.title = element_text(size=18),
        legend.key.height = unit(0.8, 'cm'),
        legend.text = element_text(size=16))
p4_leg_leg <- cowplot::get_legend(p4_leg+
                                    theme(legend.justification = 
                                            c(0,2.2)))

legend <- cowplot::plot_grid(p1_leg_leg,
                             p3_leg_leg,
                             p2_leg_leg,
                             p4_leg_leg,
                             nrow = 4)

p1 <- cowplot::plot_grid(pp1,pp2,nrow = 2,
                         rel_heights = c(4.7,5.3),
                         labels = c("(a)","(b)"),
                         vjust= 1.09,
                         label_size = 24, 
                         align = "v")
Figure_1 <- cowplot::plot_grid(p1,
                               legend,nrow = 1,
                               ncol = 2,
                               rel_widths = c(7,1.8))
Figure_1
#___
# Get the number of unique species in the 'df_tbp' data frame
Nspc <- length(unique(df_tbp$Species))

#setNames(df_clr$col.f_spc, df_clr$spcNm)
# get the unique species names from the 'df_tbp' data frame
utbpspc <- unique(df_tbp$Species)
# get the matching colors from the 'df_clr' data frame
col.f_spc.f.leg <- df_clr$col.f_spc[match(utbpspc,df_clr$spcNm)]
# Create a data frame for the legend
legend_data <- data.frame(
  Species = unique(df_tbp$shrtNm),
  color = col.f_spc.f.leg
)

# Create a ggplot object for the legend
p1_leg <- ggplot() +
  geom_line(data = legend_data, aes(x = 1, y = 1, 
                                    color = Species), size = 1, linetype = 1) +
  scale_color_manual(
    name = 'Assay',
    #values = setNames(clr04, unique(df_tbp$Species))
    values = setNames(df_clr$col.f_spc, df_clr$shrtNm)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(size = 18),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = 16, face = "italic")
  )

# Extract the legend
p1_leg_leg <- cowplot::get_legend(
  p1_leg + theme(legend.justification = c(0, 0))
)

# Create a ggplot object for the line type legend
p3_leg <- ggplot() +
  geom_line(aes(x = NA, y = NA, linetype = "ddPCR sensitivity"), linewidth = 1, color = "black") +
  geom_line(aes(x = NA, y = NA, linetype = "qPCR sensitivity"), linewidth = 1, color = "black") +
  theme_void() +
  scale_linetype_manual(
    name = "Modelled sensitivity",
    values = c("ddPCR sensitivity" = 2, "qPCR sensitivity" = 1)
  ) +
  theme(
    legend.position = "right",
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(size = 18),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = 16)
  )

# Extract the legend
p3_leg_leg <- cowplot::get_legend(
  p3_leg + theme(legend.justification = c(0, 0.7))  # Adjust vertical justification to move up
)

# Create a ggplot object for the sensitivity difference legend
p2_leg <- ggplot() +
  geom_line(aes(x = NA, y = NA, linetype = "Sensitivity difference\n(ddPCR - qPCR)"),
            size = 1, color = "black") +
  theme_void() +
  scale_linetype_manual(
    name = "",
    values = c("Sensitivity difference\n(ddPCR - qPCR)" = 6)
  ) +
  theme(
    legend.position = "right",
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(size = 1),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = 16)
  )

# Extract the legend
p2_leg_leg <- cowplot::get_legend(
  p2_leg + theme(legend.justification = c(0, 0.7))  # Adjust vertical justification to move up
)

# Create a ggplot object for the point type legend
p4_leg <- ggplot() +
  geom_point(aes(x = NA, y = NA, shape = "ddPCR standard samples"), size = 3, color = "black") +
  geom_point(aes(x = NA, y = NA, shape = "qPCR standard samples"), size = 3, color = "black") +
  theme_void() +
  scale_shape_manual(
    name = "Measured sensitivity",
    values = c("ddPCR standard samples" = 19, "qPCR standard samples" = 15)
  ) +
  theme(
    legend.position = "right",
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(size = 18),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = 16)
  )

# Extract the legend
p4_leg_leg <- cowplot::get_legend(
  p4_leg + theme(legend.justification = c(0, 0.7))  # Adjust vertical justification to move up
)

# Combine all legends with adjusted vertical positioning
legend <- cowplot::plot_grid(
  p1_leg_leg,
  p3_leg_leg,
  p2_leg_leg,
  p4_leg_leg,
  nrow = 4,
  rel_heights = c(2, 0.8, 0.8, 0.8)  # Adjust the relative heights to bring legends closer together
)

# Adjust the main plots to remove headers and ensure uniform subplot sizes
pp1 <- pp1 +
  theme(
    strip.text = element_blank(),  # Remove the headers above each subplot
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12)  # Adjust the font size of the x-axis tick labels
  )

pp2 <- pp2 +
  theme(
    strip.text = element_blank(),  # Remove the headers above each subplot
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12)  # Adjust the font size of the x-axis tick labels
  )

# Combine the main plot and the legend
p1 <- cowplot::plot_grid(pp1, pp2, nrow = 2, 
                         rel_heights = c(1, 1), 
                         labels = c("(a)", "(b)"),
                         vjust= 1.09,
                         label_size = 24, align = "v")
Figure_1 <- cowplot::plot_grid(p1, legend, nrow = 1, ncol = 2, rel_widths = c(7, 1.8))

# Display the final figure
Figure_1


# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = Figure_1, 
         # define the output filename by pasting together 
         filename = paste0(outdir15,
                           "/Fig17_Nominal_DNA_",
                           "concentration_and_detected_for_",
                           ".png"),
         width=210*2.4,height=297*1.2,
         #width=297,height=210,
         units="mm",dpi=300)
}
#____________________________________________________________________________#
#____________________________________________________________________________#

diff <- df_dqd
diff$Species <- diff$shrtNm
df_clr$Species <- df_clr$shrtNm
pp1 <-
  diff %>%
  ggplot()+
  geom_errorbar(aes(y=10^C_d,
                    x=10^C_q,
                    ymin=10^d2,
                    ymax=10^d98),
                color="grey")+
  geom_errorbar(aes(y=10^C_d,
                    x=10^C_q,
                    xmin=10^q2,
                    xmax=10^q98),
                color="grey")+
  geom_point(aes(y=10^C_d,
                 x=10^C_q,
                 color=Species))+
  scale_color_manual(
    # use the color range for the number
    # of species in the data frame
    # as inferred above
    values = setNames(df_clr$col.f_spc, df_clr$Species))+
  # c("tomato2", "deepskyblue2","orange2"))+
  geom_abline(intercept = 0,slope=1,lty=2)+
  theme_bw()+
  scale_y_log10(labels=scientific_10,
                lim=c(1e-4,1e4),
                breaks=10^c(-4,-2,0,2,4))+
  scale_x_log10(labels=scientific_10,
                lim=c(1e-5,2e2),
                breaks=10^c(-4,-2,0,2))+
  facet_wrap(~Species, nrow=3)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 19),
        plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
        axis.title.x = element_text(size = 19),
        axis.text.x = element_text(size=18),
        axis.text.y=element_text(size=14))+
  ylab("ddPCR modelled concentration \n(copies/µL)")+
  xlab("qPCR modelled concentration (copies/µL)")

#pp1


# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = pp1, 
         # define the output filenmae by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir15,
                           "/Fig18_01_ddPCR_and_qPCR_modelled_concentration",
                           ".png"),
         width=210*1.6,height=297*0.6,
         #width=297,height=210,
         units="mm",dpi=300)
}
# count the number of 'C_d' values per Species that are greater than 0
diff$Species <- diff$Species
pnts.psp <- diff %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(n = sum(C_d>0)) %>%
  dplyr::arrange(desc(n))
# get the Species name from the 'pnts.psp' data frame where the 'n'
# value is less than 3
# these needs to be subsetted from the data frame, as it is not possible
# to do a 'geom_smooth' line in ggplot2 if there are less than 3 points
spc.to.sbst <- pnts.psp %>%
  dplyr::filter(n<3) %>%
  dplyr::pull(Species)
# get the number of species in the 'diff' data frame
Nspc <- length(unique(diff$Species))

# subset the diff data frame to only include the Species in spc.to.sbst
# to produce a data frame that contains the Species that have less than 3 points
diff.sbslt3 <- diff %>%
  dplyr::filter(Species %in% spc.to.sbst)
# subset the diff data frame to include all other Species than the species in spc.to.sbst
# to produce a data frame that contains the Species that have at least 3 points
diff.sbsmt3 <- diff %>%
  dplyr::filter(!Species %in% spc.to.sbst)
# define bl3MakePlot and bm3MakePlot variables to use later on
# to determine whether to make a plot or not, setting these to FALSE
# make the default setting not to make a plot
bl3MakePlot <- FALSE
bm3MakePlot <- FALSE
# evaluate on these two data frames, whether to make a plot or not
# Check if the 'diff.sbsmt3' is empty
if(nrow(diff.sbsmt3)!=0){
  # if it is not empty, then set the 'bm3MakePlot' to TRUE
  bm3MakePlot <- TRUE}
# also evaluate on the 'diff.sbslt3' data frame
if(nrow(diff.sbslt3)!=0){
  # if it is not empty, then set the 'bl3MakePlot' to TRUE
  bl3MakePlot <- TRUE
}
# The 'bl3MakePlot' indicates whether a plot should be created for
# the Species that have less than 3 points, and 'bm3MakePlot' indicates
# whether a plot should be created for the Species that have more then 3 points or more
bl3MakePlot
bm3MakePlot
# remove any previous versions of the plots
rm(pp2.l3)
rm(pp2.m3)
# Calculate the minimum and maximum values of 10^(d98-d2)
y_values <- 10^(diff.sbsmt3$d98 - diff.sbsmt3$d2)
y_min <- min(y_values)
y_max <- max(y_values)

# Calculate the nearest powers of 10
y_min_log10 <- floor(log10(y_min))
y_max_log10 <- ceiling(log10(y_max))

# Set the limits to the nearest powers of 10
y_min_limit <- 10^y_min_log10
y_max_limit <- 10^y_max_log10

# Print the calculated min and max limits for verification
print(paste("y_min_limit:", y_min_limit))
print(paste("y_max_limit:", y_max_limit))

# if the 'bl3MakePlot' is TRUE, then prepare the pp2 plot without the geom_smooth line
# but use the 'diff.sbslt3' data frame, that contains the Species that have 
# less than 3 points
if (bl3MakePlot==T){
  pp2.l3 <-
    diff.sbslt3 %>% 
    dplyr::filter(C_d>-2) %>%
    dplyr::filter(C_q>-2) %>%
    ggplot()+
    facet_wrap(~Species, nrow=3)+
    # geom_smooth(aes(y=10^(d98-d2),x=10^C_d,color=Species),
    #             lty=2,se=T)+
    # geom_smooth(aes(y=10^(q98-q2),x=10^C_q,color=Species),
    #             lty=1,se=T)+
    geom_point(aes(y=10^(d98-d2),x=10^C_d,color=Species),
               size=2)+
    geom_point(aes(y=10^(q98-q2),x=10^C_q,color=Species),
               size=2,shape=15)+
    scale_color_manual(
      # use the color range for the number
      # of species in the data frame
      # as inferred above
      values = setNames(df_clr$col.f_spc, df_clr$Species))+
    # c("tomato2", "deepskyblue2","orange2"))+
    scale_y_log10(labels = scientific_10, limits = c(y_min_limit, y_max_limit)) +
    scale_x_log10(labels=scientific_10,lim=c(1e-5,1e2),
                  breaks=10^c(-4,-2,0,2))+
    theme_bw()+
    #facet_grid(~Species)+
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.title.y = element_text(size = 19),
          plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
          axis.title.x = element_text(size = 19),
          axis.text.x = element_text(size=18),
          axis.text.y=element_text(size=14))+
    ylab(expression("95% Credible\n Interval range"))+
    xlab("Modelled concentration (copies/µL)")
  
  pp2.l3}

# if the # if the 'bm3MakePlot' is TRUE, then prepare the pp2 plot with the 
# geom_smooth line
# but use the 'diff.sbslt3' data frame, that contains the Species that have 
# more than 3 points
if (bm3MakePlot==T){
  pp2.m3 <-
    # Create the plot with dynamic limits
    diff.sbsmt3 %>%
    dplyr::filter(C_d > -2) %>%
    dplyr::filter(C_q > -2) %>%
    ggplot() +
    facet_wrap(~Species, nrow=3) +
    geom_point(aes(y = 10^(d98-d2), x = 10^C_d, color = Species), size = 2) +
    geom_point(aes(y = 10^(q98-q2), x = 10^C_q, color = Species), size = 2, shape = 15) +
    geom_smooth(aes(y = 10^(d98-d2), x = 10^C_d, color = Species), lty = 2, se = TRUE) +
    geom_smooth(aes(y = 10^(q98-q2), x = 10^C_q, color = Species), lty = 1, se = TRUE) +
    scale_color_manual(values = setNames(df_clr$col.f_spc, df_clr$spcNm)) +
    scale_y_log10(labels = scientific_10, limits = c(y_min_limit, y_max_limit)) +
    scale_x_log10(labels = scientific_10, limits = c(1e-5, 1e2), breaks = 10^c(-4, -2, 0, 2)) +
    theme_bw() +
    #facet_grid(~Species)+
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.title.y = element_text(size = 19),
          plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
          axis.title.x = element_text(size = 19),
          axis.text.x = element_text(size=18),
          axis.text.y=element_text(size=14))+
    ylab(expression("95% Credible \nInterval range"))+
    xlab("Modelled concentration (copies/µL)")
  
  pp2.m3}
# check if the pp2.l3 and pp2.m3 exists
exists("pp2.m3")
exists("pp2.l3")

bSaveFigures.m3 <- F
bSaveFigures.l3 <- F
# evaluate on the bSaveFigures.m3 and bSaveFigures.l3 variable to determine 
# whether to save the figure or not
if (exists("pp2.l3")==T){
  bSaveFigures.l3 <- T
}
# save the figure if the above 'bSaveFigures.l3' is TRUE
if(bSaveFigures.l3==T){
  ggsave(plot = pp2.l3, 
         # define the output filename by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir15,
                           "/Fig18_02_less_than_3_ddPCR_and_qPCR_modelled_concentration",
                           ".png"),
         width=210*1.6,height=297*0.6,
         #width=297,height=210,
         units="mm",dpi=300)
}
# evaluate on the 'pp2.m3' to set the 'bSaveFigures.m3' variable
if (exists("pp2.m3")==T){
  bSaveFigures.m3 <- T
}
# save the figure if the above 'bSaveFigures.m3' is TRUE
if(bSaveFigures.m3==T){
  ggsave(plot = pp2.m3, 
         # define the output filename by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir15,
                           "/Fig18_02_more_than_3_ddPCR_and_qPCR_modelled_concentration",
                           ".png"),
         width=210*1.6,height=297*0.6,
         #width=297,height=210,
         units="mm",dpi=300)
}
# use either the pp2.m3 or the pp2.l3 to make it up for the pp2 plot
# if the 'pp2.m3' exists, then use it, otherwise use the 'pp2.l3'
if (exists("pp2.m3")==T){
  pp2 <- pp2.m3
} else {
  pp2 <- pp2.l3
}

p1_leg <- ggplot() +
  geom_line(data = legend_data, 
            aes(x = 1, y = 1, color = Species), size = 1, linetype = 1) +
  scale_color_manual(
    name = 'Assay',
    values = setNames(df_clr$col.f_spc, df_clr$Species) ) +
  
  theme_bw()+
  theme(plot.margin = margin(1, 4, 1, 0, "cm"),
        legend.key.width = unit(1.2,"cm"),
        legend.title = element_text(size=18),
        legend.key.height = unit(0.8, 'cm'),
        legend.text = element_text(size=16,face = "italic"))

p1_leg_leg <- cowplot::get_legend(p1_leg+
                                    theme(legend.justification = 
                                            c(0,0.1)))

p2_leg <- ggplot() +
  geom_line(aes(x = NA, y = NA, color = "ddPCR"), size = 1) +
  geom_line(aes(x = NA, y = NA, color = "qPCR"), size = 1) +
  theme_bw()+
  scale_color_manual(
    name="Quantification precision \n(expressed as variance)",
    breaks=c("ddPCR","qPCR"),
    values = c("black","black"),
    guide = guide_legend(
      override.aes = list(lty = c(2,1),
                          size = c(1,2)))) +
  theme(plot.margin = margin(1, 4, 1, 0, "cm"),
        legend.key.width = unit(1.2,"cm"),
        legend.title = element_text(size=18),
        legend.key.height = unit(0.8, 'cm'),
        legend.text = element_text(size=16,face = "italic"))
p2_leg_leg <- cowplot::get_legend(p2_leg+
                                    theme(legend.justification = 
                                            c(0,1.0)))

legend <- cowplot::plot_grid(p1_leg_leg,p2_leg_leg,nrow = 2)
p1 <- cowplot::plot_grid(pp1,pp2,nrow=2,
                         align = "v",
                         rel_heights = c(4,3),
                         label_size = 24,
                         labels = c("(a)","(b)"),
                         vjust= 1.2,
                         label_x = -0.003)
Figure_2 <- cowplot::plot_grid(p1,
                               legend,
                               ncol=2,
                               rel_widths = c(4,1.0))
Figure_2

# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = Figure_2, 
         # define the output filenmae by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir15,
                           "/Fig18_03_ddPCR_and_qPCR_modelled_concentration",
                           ".png"),
         width=210*2.4,height=297*0.8,
         #width=297,height=210,
         units="mm",dpi=300)
}
#___

pp1 <-
  diff %>%
  ggplot()+
  geom_errorbar(aes(y=10^C_d,
                    x=10^C_q,
                    ymin=10^d2,
                    ymax=10^d98),
                color="grey")+
  geom_errorbar(aes(y=10^C_d,
                    x=10^C_q,
                    xmin=10^q2,
                    xmax=10^q98),
                color="grey")+
  geom_point(aes(y=10^C_d,
                 x=10^C_q,
                 color=Species))+
  scale_color_manual(
    # use the color range for the number
    # of species in the data frame
    # as inferred above
    values = setNames(df_clr$col.f_spc, df_clr$Species))+
  # c("tomato2", "deepskyblue2","orange2"))+
  geom_abline(intercept = 0,slope=1,lty=2)+
  theme_bw()+
  scale_y_log10(labels=scientific_10,
                lim=c(1e-4,1e4),
                breaks=10^c(-4,-2,0,2,4))+
  scale_x_log10(labels=scientific_10,
                lim=c(1e-5,2e2),
                breaks=10^c(-4,-2,0,2))+
  facet_wrap(~Species, nrow=3)+
  theme(#legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_text(size = 17),
    plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
    axis.title.x = element_text(size = 17),
    axis.text.x = element_text(size=14),
    axis.text.y=element_text(size=14),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12, face = "italic") )+
  guides(color=guide_legend("Assay")) +
  ylab("ddPCR modelled concentration \n(copies/µL)")+
  xlab("qPCR modelled concentration (copies/µL)")

#pp1


# save the figure if the above 'bSaveFigures' is TRUE
if(bSaveFigures==T){
  ggsave(plot = pp1, 
         # define the output filenmae by pasting together 
         # qpcrrunno and qpcrrundate
         filename = paste0(outdir15,
                           "/Fig18_04_ddPCR_and_qPCR_modelled_concentration",
                           ".png"),
         width=210*1.6,height=297*0.6,
         #width=297,height=210,
         units="mm",dpi=300)
}

