library(tidyverse)
library(scales)

# --- Constants ---
genome_length <- 9719
x_axis_max <- 10000

# --- Regions ---
pol_region <- tibble(
  start = 2085,
  end = 5096
)

pol_subregions <- tibble(
  gene = c("PR", "RT", "p15", "INT"),
  start = c(2253, 2550, 3870, 4230),
  end = c(2550, 3870, 4230, 5096)
)

primer_colors <- c(
  "PR-576" = "#0D4F8B",    # Dark blue
  "RT-718" = "#33A1C9",    # Blue-cyan
  "INT-948" = "#FFCC11"    # Yellow
)

primers <- read_csv("data/regions_igs_no_gag.csv") %>%
  separate(product_coords_incl_primers, into = c("prod_start", "prod_end"), sep = "-", convert = TRUE) %>%
  mutate(pair = str_trim(pair)) %>%
  filter(pair %in% names(primer_colors))

# --- Separator lines including one before PR (at 2253) ---
separator_lines <- bind_rows(
  tibble(start = 2253, ymin = 0, ymax = 0.035),  # Dashed line before PR
  pol_subregions %>%
    filter(start != 2253) %>%  # Other subregions except PR
    transmute(start = start, ymin = 0, ymax = 0.035)
)

x_ticks <- seq(0, x_axis_max, by = 1000)
x_tick_labels <- tibble(
  x = x_ticks,
  y = -0.165,
  label = as.character(x_ticks)
)

# --- Legend data
legend_data <- tibble(
  pair = c("PR-576", "RT-718", "INT-948"),
  label = c("2253–2550", "2550–3870", "4230–5096"),  # Correct coordinates
  xstart = 7550,
  xend = 7700,
  y = 0.035 - (0.015 * (0:2)),
  color = c(primer_colors["PR-576"], primer_colors["RT-718"], primer_colors["INT-948"])
)


products_plot <- primers %>%
  mutate(y = case_when(
    pair == "PR-576" ~ -0.085,   # above genome line
    pair == "RT-718" ~ -0.113,  # below genome line but above x-axis
    pair == "INT-948" ~ -0.085   # above x-axis, same level as PR
  ))

# --- Plot ---
p <- ggplot() +
  # Pol box with blue outline (same vertical size as before)
  geom_rect(aes(xmin = pol_region$start, xmax = pol_region$end, ymin = 0, ymax = 0.035),
            fill = NA, color = primer_colors["PR-576"], size = 1) +
  
  # Smaller gene boxes inside pol, no outline, just text
  geom_rect(data = pol_subregions,
            aes(xmin = start, xmax = end, ymin = 0.005, ymax = 0.030),
            fill = NA, color = NA) +
  
  geom_text(data = pol_subregions,
            aes(x = (start + end) / 2, y = 0.0175, label = gene),
            size = 3.5, fontface = "plain", color = "black") +
  
  # Dashed vertical separator lines just outside gene boxes (including before PR)
  geom_segment(data = separator_lines,
               aes(x = start, xend = start, y = ymin, yend = ymax),
               linetype = "dashed", color = primer_colors["PR-576"], size = 0.5) +
  
  # PCR product lines colored by primer pair (normal thickness)
  geom_segment(data = products_plot,
               aes(x = prod_start, xend = prod_end, y = y, yend = y, color = pair),
               size = 1.2) +
  
  # Mapping baseline (genome line) above x-axis line
  geom_segment(aes(x = 0, xend = genome_length, y = -0.10, yend = -0.10),
               color = "gray60", size = 1.5) +
  
  # Manual x-axis line below mapping baseline
  geom_segment(aes(x = 0, xend = x_axis_max, y = -0.14, yend = -0.14),
               color = "black", size = 0.6) +
  
  # X-axis ticks every 1000
  geom_segment(data = x_tick_labels,
               aes(x = x, xend = x, y = -0.14, yend = -0.1435),
               color = "black", size = 0.4) +
  
  geom_text(data = x_tick_labels,
            aes(x = x, y = y, label = label),
            size = 4, vjust = 1.2, color = "black") +
  
  # Legend lines (thicker and longer) - no overlap with text
  geom_segment(data = legend_data,
               aes(x = xstart, xend = xend, y = y, yend = y, color = pair),
               size = 1.5,
               lineend = "round",
               show.legend = FALSE) +
  
  geom_text(data = legend_data,
            aes(x = xend + 250, y = y, label = label),
            hjust = 0,
            size = 4,
            color = "black") +
  

  # Scales & theme
  scale_x_continuous(
    limits = c(0, x_axis_max + 150),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(limits = c(-0.18, 0.06), expand = c(0, 0)) +
  
  scale_color_manual(name = "PCR products", values = primer_colors) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

print(p)


# Optionally save:

ggsave("images/img_png/pol_primer_product_mapping.png", 
       plot = p, width = 15, height = 3.2, dpi = 500)

ggsave("images/img_svg/pol_primer_product_mapping.svg", 
       plot = p, device = "svg", width = 15, height = 3.5)