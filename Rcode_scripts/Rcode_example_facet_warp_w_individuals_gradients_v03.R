library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(patchwork)
library(cowplot)  # For legend extraction

# Prepare data
di <- diamonds
di$colr <- di$color
di$sample <- paste0(di$color, "_", di$cut)
di %>% distinct(colr, sample) %>% arrange(colr, sample)

# Define the Set1 color palette for D-J colors
unique_colors <- sort(unique(di$colr))  # D to J
palette_colors <- brewer.pal(n = length(unique_colors), name = "Set1")
names(palette_colors) <- unique_colors

# Filter out zero values for log10 transformation
positive_x <- di %>% filter(x > 0)  # Only include positive values of x

# Compute the adjusted minimum and maximum for the log10 scale
min_log_x <- min(log10(positive_x$x))  # Min of log10(x) for positive x
max_log_x <- max(log10(di$x))  # Max log10(x)

# Function to create each individual plot
make_facet_plot <- function(color_label, color_hex, show_x_axis = FALSE) {
  subset_di <- di %>% filter(colr == color_label)
  
  ggplot(data = subset_di, aes(clarity, sample, fill = log10(x))) +
    geom_tile(colour = "white") +
    scale_fill_gradientn(
      colors = c("white", color_hex, "black"),
      name = "log10(x)",
      na.value = "grey50"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      legend.position = "none",  # No legend in individual subplots
      panel.spacing = unit(0.02, "lines"),
      plot.margin = margin(0, 0, 0, 0),  # No extra margin around plots
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
print(combined_plot)

# Create a dummy plot to extract the legend (with the same fill gradient)
dummy_plot <- ggplot(data = di, aes(clarity, sample, fill = log10(x))) +
  geom_tile(colour = "white") +
  scale_fill_gradientn(
    colors = c("white", "grey", "black"),
    limits = c(min_log_x, 1.5),  # Use the adjusted limits for the color scale
    breaks = c(0, 0.5, 1, 1.5),  # Log scale breaks at log10(1) = 0, log10(3) = 0.5, log10(10) = 1, log10(32) = 1.5
    labels = c("1", "3", "10", "32"),  # Actual values corresponding to the breaks
    name = "log10(x)",
    na.value = "grey50"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# Extract the legend from the dummy plot
legend <- cowplot::get_legend(dummy_plot)

# Combine the individual plots with the legend
final_plot <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(0.85, 0.15))  # Adjust width ratio

# Display final plot
final_plot

#______
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(patchwork)
library(cowplot)

# Prepare data
di <- diamonds
di$colr <- di$color
di$sample <- paste0(di$color, "_", di$cut)

# Define the Set1 color palette for D-J colors
unique_colors <- sort(unique(di$colr))  # D to J
palette_colors <- brewer.pal(n = length(unique_colors), name = "Set1")
names(palette_colors) <- unique_colors

# Filter out zero values for log10 transformation
positive_x <- di %>% filter(x > 0)  # Only include positive values of x

# Compute the adjusted minimum and maximum for the log10 scale
min_log_x <- min(log10(positive_x$x))  # Min of log10(x) for positive x
max_log_x <- max(log10(di$x))  # Max log10(x)

# Function to create each individual plot
make_facet_plot <- function(color_label, color_hex, show_x_axis = FALSE) {
  subset_di <- di %>% filter(colr == color_label)
  
  ggplot(data = subset_di, aes(clarity, sample, fill = log10(x))) +
    geom_tile(colour = "white") +
    scale_fill_gradientn(
      colors = c("white", color_hex, "black"),
      name = "log10(x)",
      na.value = "grey50"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      legend.position = "none",  # No legend in individual subplots
      panel.spacing = unit(0.02, "lines"),
      plot.margin = margin(0, 0, 0, 0),  # No extra margin around plots
      axis.text.x = if (show_x_axis) {
        element_text(angle = 90, vjust = 0.5, hjust = 1)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      strip.text = element_text(size = 14, face = "bold")  # Customize the strip labels (D, E, F, etc.)
    ) +
    scale_y_discrete(limits = rev) +
    facet_wrap(~colr, scales = "free_y", ncol = 1, strip.position = "left")  # Use facet_wrap with customized strip text
}

# Create all plots, enabling x-axis only on the last (J)
plot_list <- lapply(unique_colors, function(col) {
  make_facet_plot(col, palette_colors[col], show_x_axis = col == "J")
})

# Combine using patchwork with equal height for each subplot
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 1) &
  theme(plot.margin = margin(0, 0, 0, 0))  # Ensure no extra margin globally

# Create a dummy plot to extract the legend (with the same fill gradient)
dummy_plot <- ggplot(data = di, aes(clarity, sample, fill = log10(x))) +
  geom_tile(colour = "white") +
  scale_fill_gradientn(
    colors = c("white", "grey", "black"),
    limits = c(min_log_x, 1.5),  # Use the adjusted limits for the color scale
    breaks = c(0, 0.5, 1, 1.5),  # Log scale breaks at log10(1) = 0, log10(3) = 0.5, log10(10) = 1, log10(32) = 1.5
    labels = c("1", "3", "10", "32"),  # Actual values corresponding to the breaks
    name = "log10(x)",
    na.value = "grey50"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# Extract the legend from the dummy plot
legend <- cowplot::get_legend(dummy_plot)

# Combine the individual plots with the legend
final_plot <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(0.85, 0.15))  # Adjust width ratio

# Display final plot
final_plot
#_____

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(patchwork)
library(cowplot)

# Prepare data
di <- diamonds
di$colr <- di$color
di$sample <- paste0(di$color, "_", di$cut)

# Define the Set1 color palette for D-J colors
unique_colors <- sort(unique(di$colr))  # D to J
palette_colors <- brewer.pal(n = length(unique_colors), name = "Set1")
names(palette_colors) <- unique_colors

# Filter out zero values for log10 transformation
positive_x <- di %>% filter(x > 0)  # Only include positive values of x

# Compute the adjusted minimum and maximum for the log10 scale
min_log_x <- min(log10(positive_x$x))  # Min of log10(x) for positive x
max_log_x <- max(log10(di$x))  # Max log10(x)

# Function to create each individual plot
make_facet_plot <- function(color_label, color_hex, show_x_axis = FALSE) {
  subset_di <- di %>% filter(colr == color_label)
  
  ggplot(data = subset_di, aes(clarity, sample, fill = log10(x))) +
    geom_tile(colour = "white") +
    scale_fill_gradientn(
      colors = c("white", color_hex, "black"),
      name = "log10(x)",
      na.value = "grey50"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      legend.position = "none",  # No legend in individual subplots
      panel.spacing = unit(0.02, "lines"),
      plot.margin = margin(0, 0, 0, 0),  # No extra margin around plots
      axis.text.x = if (show_x_axis) {
        element_text(angle = 0, vjust = 0.5, hjust = 1)  # Rotate x-axis labels by 90 degrees
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      axis.text.y = element_blank(),  # Remove the y-axis text (sample labels)
      axis.ticks.y = element_blank(),  # Optional: Remove the y-axis ticks
      strip.text.y.left = element_text(size = 14, face = "bold", angle = 0, vjust = 0.5),  # Rotate the facet labels for y-axis
      strip.text.x = element_text(size = 14, face = "bold", angle = 90, vjust = 0.5)  # Optionally rotate facet labels on x-axis (for clarity)
    ) +
    scale_y_discrete(limits = rev) +
    facet_wrap(~colr, scales = "free_y", ncol = 1, strip.position = "left")  # Use facet_wrap with customized strip text
}

# Create all plots, enabling x-axis only on the last (J)
plot_list <- lapply(unique_colors, function(col) {
  make_facet_plot(col, palette_colors[col], show_x_axis = col == "J")
})

# Combine using patchwork with equal height for each subplot
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 1) &
  theme(plot.margin = margin(0, 0, 0, 0))  # Ensure no extra margin globally

# Create a dummy plot to extract the legend (with the same fill gradient)
dummy_plot <- ggplot(data = di, aes(clarity, sample, fill = log10(x))) +
  geom_tile(colour = "white") +
  scale_fill_gradientn(
    colors = c("white", "grey", "black"),
    limits = c(min_log_x, 1.5),  # Use the adjusted limits for the color scale
    breaks = c(0, 0.5, 1, 1.5),  # Log scale breaks at log10(1) = 0, log10(3) = 0.5, log10(10) = 1, log10(32) = 1.5
    labels = c("1", "3", "10", "32"),  # Actual values corresponding to the breaks
    name = "log10(x)",
    na.value = "grey50"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# Extract the legend from the dummy plot
legend <- cowplot::get_legend(dummy_plot)

# Combine the individual plots with the legend
final_plot <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(0.85, 0.15))  # Adjust width ratio

# Display final plot
final_plot

