library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(patchwork)

# Prepare data
di <- diamonds
di$colr <- di$color
di$sample <- paste0(di$color, "_", di$cut)

# Define the Set1 color palette for D-J colors
unique_colors <- sort(unique(di$colr))  # D to J
palette_colors <- brewer.pal(n = length(unique_colors), name = "Set1")
names(palette_colors) <- unique_colors
max(di$x)
min(log10(di$x))
min(di$x)
# Function to create each individual plot
make_facet_plot <- function(color_label, color_hex, show_x_axis = FALSE) {
  subset_di <- di %>% filter(colr == color_label)
  
  ggplot(data = subset_di, aes(clarity, sample, fill = log10(x))) +
    geom_tile(colour = "white") +
    scale_fill_gradientn(colors = c("white", color_hex, "black"),
                         name = "log10(x)",
                         na.value = "grey50") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.02, "lines"),
      plot.margin = margin(0, 0, 0, 0),  # Remove margin around each subplot
      axis.text.x = if (show_x_axis) {
        element_text(angle = 90, vjust = 0.5, hjust = 1)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) element_line() else element_blank()
    ) +
    scale_y_discrete(limits = rev)
}

# Create all plots, enabling x-axis only on the last (J)
plot_list <- lapply(unique_colors, function(col) {
  make_facet_plot(col, palette_colors[col], show_x_axis = col == "J")
})

# Combine using patchwork with equal height for each subplot
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 1) &
  theme(plot.margin = margin(0, 0, 0, 0))  # Ensure no extra margin globally

combined_plot
