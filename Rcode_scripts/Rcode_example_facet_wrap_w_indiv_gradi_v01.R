di<- diamonds
di$colr<- di$color
length(unique(di$color))
di$sample <- paste0(di$color,"_", di$cut)
# make a color palette of colors that can reflect the individual colors in 'di$colr'
library(ggplot2)
library(RColorBrewer)
# Define a color palette
m_plt <- brewer.pal(n = 7, name = "Set1")
# for the distinct elements in 'di$colr', get a color from the 'm_plt' palette
di <- di %>%
  mutate(colr = factor(colr, 
                       levels = unique(colr), labels = m_plt[1:length(unique(colr))]))
# Create the heatmap
gplot4 <- ggplot(data = di, 
                 aes(clarity, sample, 
                     fill = log10(x))) +
  geom_tile(colour = "white") +
  facet_wrap(~colr, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(strip.text = element_blank()) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradientn(colors = c("white", m_plt, "black"),
                       name = "log10(x)", 
                       na.value = "grey50") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position="none") +
  theme(panel.spacing = unit(0.02, "lines"))
gplot4