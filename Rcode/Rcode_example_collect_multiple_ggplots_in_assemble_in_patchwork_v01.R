#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#dev.off()

########################################################################################
#remove everything in the working environment, without a warning!!
rm(list=ls())


# example on collection of ggplots
# while iterating over subsetted  data frames 
# i.e. -  iterating over different input files
# for getting different  data frames 
library(ggplot2)
mtcars2 <- mtcars
mtcars2 <- mtcars2[grepl("Merc",row.names(mtcars2)),]
mtcars2$wellID <- mtcars2$setNm
mtcars2$carNm <- row.names(mtcars2)
mtcars2 <- mtcars2[sample(nrow(mtcars2), 300, replace = T), ]
#
mtcars2$v0 <- rnorm(nrow(mtcars2), mean = 20, sd = 1)*mtcars2$disp
mtcars2$v1 <- rnorm(nrow(mtcars2), mean = 30, sd = 1)*mtcars2$disp
mtcars2$v2 <- rnorm(nrow(mtcars2), mean = 40, sd = 1)*mtcars2$disp
mtcars2$Amplitudes <- mtcars2$v2
mtcars2$setNm<- sample((LETTERS[1:4]), 300, replace = T)
mtcars2$wellID <- mtcars2$carNm
mtcars2$stdllvl <- mtcars2$carNm
ch <- 1
chCol <- c("dodgerblue3", "forestgreen")
channel <- ch
# make a vector of the species abbreviations to iterate over
set_of_spc.Abbr <- c("setA","setB","setC","setD")
# get the total number of elements in the vector containing species 
# abbreviations names
nspcAbr <- length(set_of_spc.Abbr)
#nspcAbr <- 1
# make an empty list that plots can be collected in
lst_plts <- list()

# iterate over elements in the vector of numbers for species abbreviations
for (spcAbbr.no in seq(1,nspcAbr,1))
{
  print(spcAbbr.no)#}
  # get the species abbreviation that matches the species number
  spc.Abbr_iso <- set_of_spc.Abbr[spcAbbr.no]
  spc.Abbr_iso <- gsub("set","",spc.Abbr_iso)
  PDCe_dd.01 <- mtcars2[(mtcars2$setNm==spc.Abbr_iso),]
plt02 <- ggplot(data = PDCe_dd.01) + 
  theme_bw() +
  geom_point(aes(x = seq_len(nrow(PDCe_dd.01)),
                 y = Amplitudes, group = wellID), 
                 color = "blue",size = 1) + 
  labs(x = "Event", y = "Amplitude", color = NULL) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none") + 
  # geom_hline(data = thrDfCh.01, 
  #            aes(yintercept = thrCh), 
  #            col = "magenta") + 
  scale_color_manual(labels = c("neg", 
                                "pos"), 
                     values = c("gray50", 
                                chCol[channel])) +
  # Add a label that can hold the dilution step information
  # https://stackoverflow.com/questions/36089962/ggplot2-change-strip-text-position-in-facet-grid-plot
  geom_text(aes(label = stdllvl),
            col="black",
            x = Inf, y = Inf, 
            hjust = 1.5, vjust = 1.5, check_overlap = TRUE) +
  # reduce the spacing between the panels in facet wrap
  theme(panel.spacing = unit(0.1, "lines")) +
  facet_wrap(~wellID,  ncol = 8) +
  ggtitle(spc.Abbr_iso)

#plt02
# place the plot in the list that collects the plots
lst_plts[[spcAbbr.no]] <- plt02
}

library(ggplot2)
library(gridExtra)
library(ggpubr)
lst_plts[[1]]
lst_plts[[4]]
ggpubr::ggarrange(plotlist = lst_plts,
                  nrow = 2,
                  ncol = ceiling(length(lst_plts)/2))