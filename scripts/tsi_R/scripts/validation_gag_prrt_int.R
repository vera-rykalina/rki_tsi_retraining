
# 1. Libraries ------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)
library(ggpmisc)
library(greekLetters)

# 2. Export dataframe -----------------------------------------------------


### GAG PRRT INT ###
hivtime_gag_prrt_int_df <- read_excel("data/validation_dataset_gag_prrt_int.xlsx",sheet = "GAG_PRRT_INT")
head(hivtime_gag_prrt_int_df) 
nrow(hivtime_gag_prrt_int_df)

View(hivtime_gag_prrt_int_df)

# 3. Filter dataframe -----------------------------------------------------

validation_filtered <- hivtime_gag_prrt_int_df |>
  select(scount, known_tsi_months, est_tsi_months, subtype_prrt) |> 
  mutate(tsi_diff = known_tsi_months - est_tsi_months) |> 
  mutate(tsi_levels = ifelse(known_tsi_months <= 12, "0-12", 
                             ifelse(known_tsi_months > 12 & known_tsi_months <= 24, "12-24", 
                                    ifelse(known_tsi_months > 24 & known_tsi_months <= 36, "24-36", 
                                           ifelse(known_tsi_months > 36 & known_tsi_months <= 48, "36-48",
                                                  "48-60")))))

# 4. Create a boxplot -----------------------------------------------------

# Sample size
(tsi_groups <- validation_filtered |>
   group_by(tsi_levels) |> 
   summarise(n = n()))

# Subtype groups
(tsi_subtypes <- validation_filtered |> 
    group_by(subtype_prrt) |> 
    summarise(n = n()))

subtype <- c("B", "non-B")
samples <- c(61, 12)
subtypes <- data.frame(Subtype = subtype, Count = samples) 

# Join the data with the group sizes and modify tsi_levels
(validation_filtered |> 
    left_join(tsi_groups, by = "tsi_levels") |> 
    mutate(tsi_levels_label = paste0(tsi_levels, "\n", "(n=", n,")")) |> 
    ggplot(aes(x = factor(tsi_levels_label, 
                          levels = c("0-12\n(n=53)", "12-24\n(n=6)",
                                     "24-36\n(n=4)", "36-48\n(n=9)",
                                     "48-60\n(n=1)")),
               y = tsi_diff)) +
    
    # Horizontal line at y = 0
    geom_hline(yintercept = 0, size = 0.6, color = "#3A5894") +
    
    # Error bars with colors corresponding to tsi_levels
    stat_boxplot(lwd = 0.7, geom = 'errorbar', 
                 aes(color = factor(tsi_levels_label))) +
    
    # Boxplots with colors corresponding to tsi_levels
    geom_boxplot(lwd = 0.7,
                 aes(color = factor(tsi_levels_label)), 
                 outlier.shape = 21,
                 outlier.size = 0.85,
                 outlier.stroke = 0.85) +
    
    # Custom colors for each tsi_levels
    scale_color_manual(values = c("#0D4F8B", "#33A1C9", "#33A1DE", 
                                  "#FFCC11", "#236B8E", "#0099CC"),
                       name = "Group size", 
                       labels = c("0-12 (n=53)", "12-24 (n=6)", "24-36 (n=4)", 
                                  "36-48 (n=9)", "48-60 (n=1)")) +
    
    # Jitter for individual data points
    geom_jitter(size = 0.5, alpha = 0.8, color = "grey15") +
    
    # Labels and axis formatting
    labs(x = 'Months', y = 'Known TSI - Estimated TSI') +
    scale_y_continuous(limits = c(-80, 80), breaks = seq(-80, 80, 20)) +
    
    # Add a text box
    geom_label(label = "Subtype: B (n=61), non-B (n=12)", 
               x = 1.5, y = 70,
               label.padding = unit(0.55, "lines"),
               label.size = 0.35,
               color = "grey10",
               fill = "white") +
    # Theme settings
    theme_minimal() +
    theme(text = element_text(size = 16, color = "grey10"), 
          panel.background = element_rect(fill = "white", colour = "grey80"), 
          axis.title = element_text(size = 16, color = "grey10"), 
          legend.text = element_text(size = 16),
          legend.position = "none"))


# Annotation of table
# annotate(geom = "table", label = list(subtypes), 
#         x = 0.7, y = 80, color = "grey10")

# Save the plot
ggsave("images/img_png/validation_hivtime_300_dpi.png", 
       units = "cm", width = 24, height = 12, dpi = 500)



# 5. Known vs Estimated (GAG PRRT INT) ---------------------------------------------

# Contingency table: known vs estimated (threshold is 12 months)
sk_cutoff_12 <- hivtime_gag_prrt_int_df |> 
  select(scount, known_tsi_months, est_tsi_months) |>
  mutate(est_tsi = ifelse(est_tsi_months <= 12,"recent", "long-term")) |> 
  mutate(known_tsi = ifelse(known_tsi_months <= 12,"recent", "long-term"))

sk_cutoff_12$est_tsi = paste("Estimated =", sk_cutoff_12$est_tsi)
sk_cutoff_12$known_tsi = paste("Known =", sk_cutoff_12$known_tsi)
table(sk_cutoff_12$est_tsi, sk_cutoff_12$known_tsi)

(n_sk = 18 + 2 + 17 + 36) # total number of samples 73

# Accuracy (%) known vs estimated - 74%
# Accuracy = (true positives + true negatives) / (total)
sk_hivtime_accuracy_12 = (18 + 36) / (18 + 36 + 2 + 17) * 100
sk_hivtime_accuracy_12

# FRR (False Positive)/(False Positives + True Negatives)
(sk_FRR_12 <- 2 / (2 + 18) * 100) # 10.0%
# TRR (True Positive)/(True Positives + False Negatives)
(sk_TRR_12 <- 36 / (36 + 17) * 100) # 67.9%


# 6  Known vs Estimated (GAG PRRT INT) ---------------------------------------------

# Contingency table: known vs estimated (threshold is 6 months)
sk_cutoff_6 <- hivtime_gag_prrt_int_df |> 
  select(scount, known_tsi_months, est_tsi_months) |>
  mutate(est_tsi = ifelse(est_tsi_months <= 6,"recent", "long-term")) |> 
  mutate(known_tsi = ifelse(known_tsi_months <= 6,"recent", "long-term"))

sk_cutoff_6$est_tsi = paste("Estimated =", sk_cutoff_6$est_tsi)
sk_cutoff_6$known_tsi = paste("Known =", sk_cutoff_6$known_tsi)
table(sk_cutoff_6$est_tsi, sk_cutoff_6$known_tsi)

(n_sk = 23 + 0 + 16 + 34) # total number of samples 73

# Accuracy (%) known vs estimated - 53.4%
# Accuracy = (true positives + true negatives) / (total)
sk_hivtime_accuracy_6 = (16 + 23) / (16 + 23 + 0 + 34) * 100
sk_hivtime_accuracy_6

# FRR (False Positive)/(False Positives + True Negatives)
(sk_FRR_6 <- 0 / (0 + 23) * 100) # 0.0%
# TRR (True Positive)/(True Positives + False Negatives)
(sk_TRR_6 <- 16 / (16 + 34) * 100) # 32.0%




# 7. Scatter plots TSI Difference -----------------------------------------

### GAG PRRT INT ###
hivtime_sk_diff <- hivtime_gag_prrt_int_df |> 
  select(scount, known_tsi_months, est_tsi_months) |> 
  mutate(tsi_diff = known_tsi_months - est_tsi_months) |>
  mutate(diff_abs = abs(tsi_diff)) |> 
  mutate(tsi_levels = ifelse(diff_abs <= 1, 
                             paste("0 <", greeks("Delta"), "\u2264", "1"), 
                             ifelse(diff_abs > 1 & diff_abs <= 3, 
                                    paste("1 <", greeks("Delta"), "\u2264", "3"), 
                                    ifelse(diff_abs > 3 & diff_abs <= 6, 
                                           paste("3 <", greeks("Delta"), "\u2264", "6"), 
                                           ifelse(diff_abs > 6 & diff_abs <= 12, 
                                                  paste("6 <", greeks("Delta"), "\u2264", "12"),
                                                  ifelse(diff_abs > 12 & diff_abs <= 24, 
                                                         paste("12 <", greeks("Delta"), "\u2264", "24"),
                                                         paste(greeks("Delta"), "> 24"))))))) 

(hivtime_sk_diff_plot <- hivtime_sk_diff %>%
    ggplot() +
    geom_point(aes(x = known_tsi_months, y = est_tsi_months, 
                   color = factor(tsi_levels, 
                                  levels = c(paste("0 <", greeks("Delta"), "\u2264", "1"),
                                             paste("1 <", greeks("Delta"), "\u2264", "3"),
                                             paste("3 <", greeks("Delta"), "\u2264", "6"), 
                                             paste("6 <", greeks("Delta"), "\u2264", "12"),
                                             paste("12 <", greeks("Delta"), "\u2264", "24"),
                                             paste(greeks("Delta"), "> 24"))))) +
    labs(x = "Known TSI", y = "Estimated TSI") +
    labs(title = expression(italic("gag, prrt, int")*" (recency threshlod: 12 months)")) +
    scale_x_continuous(limits = c(0, 120), breaks = seq(0,140,20)) +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0,140,20)) +
    scale_color_manual(name = "TSI difference", 
                       values = c("#FFB90F", "#778899", "#5D92B1", 
                                  "#104E8B", "#293352", "#6B4226")) +
    # Add a text box
    geom_label(label = "Accuracy: 74.0%, FRR: 10.0%, TRR: 67.9%", 
               x = 55, y = 110,
               label.padding = unit(0.55, "lines"),
               label.size = 0.25,
               color = "grey10",
               fill = "white") +
    theme_minimal() +
    theme(text = element_text(size = 12, color = "grey15"), 
          panel.background = element_rect(fill = "white", colour = "grey80"), 
          axis.title = element_text(size = 12, 
                                    color = "grey15"),
          legend.text = element_text(size = 10), 
          legend.background = element_rect(color = "white")) +
    theme(legend.position = c(0.70, 0.35), 
          legend.key = element_rect(color = "grey15")))

# Save the plot
ggsave("images/img_png/hivtime_sk_diff_400_dpi.png", 
       units = "cm", width = 12, height = 12, dpi = 400)



## Training dataset ###

# Load your real data
metadata_df <- read.csv("data/metadata.csv")

# Create a simplified subtype group: "B" or "non-B"
metadata_df <- metadata_df %>%
  mutate(subtype_group = ifelse(subtype_prrt == "B", "B", "non-B"))

# Create seroconverter interval (TSI) buckets
metadata_df <- metadata_df %>%
  mutate(tsi_period = cut(known_tsi_months,
                          breaks = c(0, 12, 24, 36, 48, 60, Inf),
                          labels = c("0-12", "12-24", "24-36", "36-48", "48-60", ">60"),
                          right = FALSE))

# Count samples per TSI period
count_data <- metadata_df %>%
  filter(!is.na(tsi_period)) %>%
  count(tsi_period)

# Calculate totals for label
total_samples <- nrow(metadata_df)
total_B <- sum(metadata_df$subtype_group == "B", na.rm = TRUE)
total_nonB <- sum(metadata_df$subtype_group == "non-B", na.rm = TRUE)

# Create one-line label text
label_text <- paste0(
  "Total samples: ", total_samples, "; ",
  "Subtype: B (n=", total_B, "), non-B (n=", total_nonB, "), NA (n=1)"
)

# Your color palette for bars
my_colors <- c("#FFB90F", "#778899", "#5D92B1", "#104E8B", "#293352", "#6B4226")

# Plot
p <- ggplot(count_data, aes(x = tsi_period, y = n, fill = tsi_period)) +
  geom_bar(stat = "identity", color = "grey15") +
  geom_text(aes(label = n), vjust = -0.5, size = 5, color = "grey15") +
  scale_fill_manual(values = my_colors, guide = "none") +
  labs(x = "Dutation of infection (months)", y = "Number of Samples",
       title = "Selected Samples for Retraining Dataset (full length)") +
  theme_minimal() +
  theme(text = element_text(size = 12, color = "grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        axis.title = element_text(size = 12, color = "grey15"),
        axis.text = element_text(color = "grey15"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("label", x = 4.8, y = max(count_data$n)*1.1,
           label = label_text,
           fill = "white", color = "grey15", size = 4, 
           label.padding = unit(0.4, "lines"), label.size = 0.25) +
  scale_y_continuous(limits = c(0, 120))


# Print plot
print(p)


# Save the plot
ggsave("images/img_png/training_dataset_500_dpi.png", 
       units = "cm", width = 24, height = 12, dpi = 500)




## Validation dataset ##
validation_df <- read_excel("data/validation_dataset_gag_prrt_int.xlsx",sheet = "GAG_PRRT_INT")

# Create a simplified subtype group: "B" or "non-B"
validation_df <- validation_df %>%
  mutate(subtype_group = ifelse(subtype_prrt == "B", "B", "non-B"))

# Create seroconverter interval (TSI) buckets
validation_df <- validation_df %>%
  mutate(tsi_period = cut(known_tsi_months,
                          breaks = c(0, 12, 24, 36, 48, 60),
                          labels = c("0-12", "12-24", "24-36", "36-48", "48-60"),
                          right = FALSE))

# Count samples per TSI period
validation_count_data <- validation_df %>%
  filter(!is.na(tsi_period)) %>%
  count(tsi_period)

# Calculate totals for label
validation_total_samples <- nrow(validation_df)
validation_total_B <- sum(validation_df$subtype_group == "B", na.rm = TRUE)
validtiona_total_nonB <- sum(validation_df$subtype_group == "non-B", na.rm = TRUE)

# Create one-line label text
label_text <- paste0(
  "Total samples: ", validation_total_samples, "; ",
  "Subtype: B (n=", validation_total_B, "), non-B (n=", total_nonB, ")"
)

# Your color palette for bars
my_colors <- c("#FFB90F", "#778899", "#5D92B1", "#104E8B", "#293352")

# Plot
validation_plot <- ggplot(validation_count_data, 
                          aes(x = tsi_period, y = n, fill = tsi_period)) +
  geom_bar(stat = "identity", color = "grey15") +
  geom_text(aes(label = n), vjust = -0.5, size = 4, color = "grey15") +
  scale_fill_manual(values = my_colors, guide = "none") +
  labs(x = "Dutation of infection (months)", y = "Number of Samples",
       title = "Samples for Validation (partial genomes)") +
  theme_minimal() +
  theme(text = element_text(size = 12, color = "grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        axis.title = element_text(size = 12, color = "grey15"),
        axis.text = element_text(color = "grey15", size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("label", x = 3, y = 90,
           label = label_text,
           fill = "white", color = "grey15", size = 4, 
           label.padding = unit(0.4, "lines"), label.size = 0.3) +
  scale_y_continuous(limits = c(0, 100))


# Print plot
print(validation_plot)


# Save the plot
ggsave("images/img_png/validation_dataset_500_dpi.png", 
       units = "cm", width = 12, height = 12, dpi = 500)


### MAE ###
### GAG PRRT INT ###
mae <- hivtime_gag_prrt_int_df |> 
  select(scount, known_tsi_months, est_tsi_months) |> 
  mutate(tsi_diff = known_tsi_months - est_tsi_months) |>
  mutate(diff_abs = abs(tsi_diff)) |> 
  summarise(mae = mean(diff_abs))
