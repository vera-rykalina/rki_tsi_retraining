
# 1. Libraries ------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)
#library(extrafont)
#library(maps)
#library(mapdata)
#library(devtools)
#library(httr)
#library(rkicolors)



#font_import(prompt = FALSE) 
#loadfonts()

# 2. Export dataframe -----------------------------------------------------
hivtime_df <- read_excel("data/hiv_phylotsi_validation.xlsx",sheet = "HIVtime_v2")
head(hivtime_df) 
nrow(hivtime_df)


# 3. Filter dataframe -----------------------------------------------------

hivtime_filtered <- hivtime_df |>
  rename(known_tsi = duration_of_infection_months) |> 
  rename(est_tsi = hivtime_v2_months) |> 
  select(scount, known_tsi, est_tsi, protocol) |> 
  filter(!scount %in% "14-00725") |> 
  filter(!protocol %in% c(paste0("FL", seq(1:9)))) |>
  mutate(tsi_diff = known_tsi - est_tsi) |> 
  mutate(tsi_levels = ifelse(known_tsi <= 12, "0-12", 
                      ifelse(known_tsi > 12 & known_tsi <= 24, "12-24", 
                      ifelse(known_tsi > 24 & known_tsi <= 36, "24-36", 
                      ifelse(known_tsi > 36 & known_tsi <= 48, "36-48",
                      ifelse(known_tsi > 48 & known_tsi <= 60, "48-60", 
                                                              ">60"))))))



# 4. Create a boxplot -----------------------------------------------------

(hivtime_filtered |> 
  group_by(tsi_levels) |> 
  ggplot(aes(x = factor(tsi_levels, 
                        levels = c("0-12","12-24","24-36",
                                   "36-48","48-60",">60")), 
             y = tsi_diff)) +
  geom_hline(yintercept = 0, size = 1, color = "#3A5894") +
  stat_boxplot(lwd = 2, geom = 'errorbar', 
               aes(color = factor(tsi_levels, 
                                  levels = c("0-12","12-24","24-36",
                                             "36-48","48-60",">60")))) +
  geom_boxplot(lwd = 1.5,
               aes(color = factor(tsi_levels, 
                                levels = c("0-12","12-24","24-36",
                                           "36-48","48-60",">60"))), 
                                outlier.shape = 21,
                                outlier.size = 2,
                                outlier.stroke = 1) +
  scale_color_manual(label = c("0-12","12-24","24-36",
                               "36-48","48-60",">60"), 
                     name = "TSI Time", 
                     values = c("#0D4F8B", "#33A1C9", "#33A1DE", 
                                "#FFCC11", "#236B8E", "#0099CC")) +
 
  geom_jitter(size = 1.4, alpha = 0.8, color = "grey15") +
  labs(x = 'Months', y = 'Known TSI - Estimated TSI\n(months)') +
  scale_y_continuous(limits = c(-90, 90), breaks = seq(-90,90,30)) +
  theme_minimal() +
  #theme_rki(base_family = "ScalaSansLF-Regular") +
  theme(text = element_text(size   = 30, 
                            color  = "grey10",
                            family = "ScalaSansLF-Regular"), 
        panel.background = element_rect(fill = "white", colour = "grey80"), 
        axis.title       = element_text(size = 30, color = "grey10"), 
        legend.text      = element_text(size = 30)) +
  theme(legend.position  = "below"))


ggsave("images/img_png/hivtime_boxplot_poster_700_dpi.png", 
       units = "cm", width = 36, height = 18, dpi = 700)

#fonttable() |>  
#  as.data.frame() |>  
#  filter(grepl("Scala", FullName)) |> 
#  select(FamilyName,FullName,Bold,Italic,Symbol)




# Contingency table: avidity vs estimated (threshold is 6 months)
avidity_contingency <- hivtime_df |> 
  rename(avidity = biorad_avidity) |> 
  rename(est_tsi = hivtime_v2_months) |>
  rename(known_tsi = duration_of_infection_months) |> 
  select(scount, protocol, known_tsi, est_tsi, avidity) |> 
  filter(!scount %in% "14-00725") |> 
  filter(!protocol %in% c(paste0("FL", seq(1:9)))) |>
  drop_na() |> 
  mutate(avidity = ifelse(avidity <= 70,"recent", "long-term")) |> 
  mutate(est_tsi = ifelse(est_tsi <= 6,"recent", "long-term")) |> 
  mutate(known_tsi = ifelse(known_tsi <= 6,"recent", "long-term"))

avidity_contingency$avidity = paste("Avidity =", avidity_contingency$avidity)
avidity_contingency$est_tsi = paste("Estimated =", avidity_contingency$est_tsi)
avidity_contingency$known_tsi = paste("Known =", avidity_contingency$known_tsi)
table(avidity_contingency$avidity, avidity_contingency$known_tsi)
table(avidity_contingency$est_tsi, avidity_contingency$known_tsi)


# Accuracy (%) known vs estimated (no NA, no FL1:FL9) - 89%
# Accuracy = (true positives + true negatives) / (total)
hivtime_accuracy = (193 + 27) / (193 + 27 + 24 + 2) * 100
hivtime_accuracy


# Accuracy (%) known vs avidity (no NA, no FL1:FL9) - 89,8%
# Accuracy = (true positives + true negatives) / (total)
avidity_accuracy = (181 + 40) / (181 + 40 + 11 + 14) * 100
avidity_accuracy
