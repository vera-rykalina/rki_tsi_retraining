library(ggplot2)
library(dplyr)
library(readr)

# Load your CSV (keep all rows first, including summary rows)
df_raw <- read_csv("data/metrics.csv", col_types = cols())

# Skip first 3 summary rows (MAE, RMSE, R2)
df_plot <- df_raw[-(1:3), ]

# Convert numeric columns safely
df_plot <- df_plot %>%
  mutate(
    true_num = as.numeric(true),
    pred_num = as.numeric(pred),
    lower_95_num = as.numeric(lower_95),
    upper_95_num = as.numeric(upper_95)
  )

# Calculate stats (for the CV metrics label box)
SS_res <- sum((df_plot$true_num - df_plot$pred_num)^2)
SS_tot <- sum((df_plot$true_num - mean(df_plot$true_num))^2)
R2 <- 1 - SS_res / SS_tot
MAE <- mean(abs(df_plot$true_num - df_plot$pred_num))
RMSE <- sqrt(mean((df_plot$true_num - df_plot$pred_num)^2))

label_text <- sprintf(
  "Linear:\nR² = %.3f\nMAE = %.3f\nRMSE = %.3f",
  R2, MAE, RMSE
)

# Compute plot limits
lims <- range(c(df_plot$true_num, df_plot$pred_num, df_plot$lower_95_num, df_plot$upper_95_num), na.rm = TRUE)
margin <- 0.05 * diff(lims)


# Plot
ggplot(df_plot, aes(x = true_num, y = pred_num)) +
  # Ribbon: 95% prediction interval
  geom_ribbon(aes(ymin = lower_95_num, ymax = upper_95_num), fill = "#33A1C9", alpha = 0.2) +
  # Points
  geom_point(color = "#0D4F8B", size = 2, alpha = 0.5, shape = 16) +
  # Identity line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#FFCC11", size = 1) +
  # CV metrics label box (top-left)
  geom_label(
    label = label_text,
    x = lims[1] + margin,
    y = lims[2] - margin,
    hjust = 0,
    vjust = 1,
    label.padding = unit(0.35, "lines"),
    label.size = 0.25,
    color = "grey10",
    fill = "white"
  ) +
  # Optional ribbon annotation (top-right)
  annotate(
    "text",
    x = lims[2] - margin, 
    y = 6.25,
    label = "95% Prediction Interval (nested CV folds)",
    color = "grey40",
    size = 3.5,
    hjust = 1,
    fontface = "italic"
  ) +
  # Fixed aspect ratio and axis limits
  coord_fixed(xlim = c(-0.5, 6.5), ylim = c(-0.5, 6.5), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 6.5, by = 1)) +
  scale_y_continuous(breaks = seq(0, 6.5, by = 1)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12, color = "grey15"),
    axis.text = element_text(color = "grey15"),
    legend.text = element_text(size = 10),
    legend.background = element_rect(color = "white"),
    plot.margin = margin(10, 10, 10, 10) # reduce extra white space
  ) +
  labs(
    title = "Known vs Estimated TSI (Nested CV)",
    x = "Known TSI (years)",
    y = "Estimated TSI (years)"
  )

# Save plot
ggsave(
  "images/img_png/true_vs_predicted_ribbon_linear.png",
  units = "cm",
  width = 12,
  height = 12,
  dpi = 400
)
