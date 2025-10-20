library(ggplot2)
library(scales)

# Sample data
set.seed(123)
df <- data.frame(
  x = 10^(runif(20, -2, 5)),
  y = 10^(runif(20, -2, 5)),
  xmin = 10^(runif(20, -2.1, 4.9)),
  xmax = 10^(runif(20, -1.9, 5.1)),
  ymin = 10^(runif(20, -2.1, 4.9)),
  ymax = 10^(runif(20, -1.9, 5.1)),
  group = factor(sample(1:2, 20, replace = TRUE))
)

# Plot
ggplot(df, aes(x = x, y = y)) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    orientation = "x",
    width = 0.1
  ) +
  geom_errorbar(
    aes(xmin = xmin, xmax = xmax),
    orientation = "y",
    width = 0.1
  ) +
  geom_point(
    aes(fill = group),
    shape = 21, size = 3, stroke = 1
  ) +
  scale_x_log10(limits = c(1e-2, 1e5)) +
  scale_y_log10(limits = c(1e-2, 1e5)) +
  annotation_logticks() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal()
