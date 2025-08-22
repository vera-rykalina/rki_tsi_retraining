
# Libraries ---------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(tidyr)
library(ggsci)
library(gg.gap)
library(RColorBrewer)
library(ggpubr)
library(psych)
library(forcats)
library(ggbreak)

# Export original table
tsi_df <- read_excel("data/hiv_phylotsi_validation.xlsx",sheet = "HIVPhyloTSI")

# Combine viral load from two columns
tsi_df$viral_load <- ifelse(is.na(tsi_df$seVLzuBlutDat), tsi_df$seVLzuDatNext, tsi_df$seVLzuBlutDat)

head(tsi_df) 
nrow(tsi_df)
range(tsi_df$duration_of_infection_months)
range(tsi_df$duration_of_infection_days)

tsi_no_na <- tsi_df |>
  rename(known_tsi = duration_of_infection_months) |> 
  rename(est_tsi = hivphylotsi_months) |> 
  select(scount, known_tsi, est_tsi, viral_load) |> 
  drop_na() |> 
  mutate(tsi_diff = known_tsi - est_tsi) |> 
  mutate(tsi_levels = ifelse(known_tsi <= 6, "0-6",  
                      ifelse(known_tsi > 6  & known_tsi  <= 12, "6-12",
                      ifelse(known_tsi > 12 & known_tsi <= 24, "12-24", 
                      ifelse(known_tsi > 24 & known_tsi <= 36, "24-36", 
                      ifelse(known_tsi > 36 & known_tsi <= 48, "36-48",
                      ifelse(known_tsi > 48 & known_tsi <= 60, "48-60", 
                                                               ">60")))))))


# Ranges
range(tsi_no_na$known_tsi) # 0.0 107.6
range(tsi_no_na$est_tsi) # 1.47 114.61
range(tsi_no_na$tsi_diff) # -85.16000  62.89667

# MAE (13.8 months)
MAE = sum(abs(tsi_no_na$tsi_diff)) / nrow(tsi_no_na)
MAE

# MAE without 3 outliers (13.2 months)
tsi_no_outlier <- tsi_no_na %>%
  filter(!scount %in% c("04-00359", "10-01035", "07-01047"))

MAE_OUT = sum(abs(tsi_no_outlier$tsi_diff))/nrow(tsi_no_outlier)
MAE_OUT

# MAE without >60 group (10.9)
tsi_no_60 <- tsi_no_na %>%
  filter(tsi_levels != ">60")

MAE_60_OUT = sum(abs(tsi_no_60$tsi_diff))/nrow(tsi_no_60)
MAE_60_OUT


# MAE up to a year (12 months, 9.9)
tsi_only_12month <- tsi_no_na %>%
  filter(tsi_levels %in% c("0-6", "6-12"))

MAE_1year = sum(abs(tsi_only_12month$tsi_diff))/nrow(tsi_only_12month)
MAE_1year


# General descriptive statistics
psych::describe(tsi_no_na)
summary(tsi_no_na)
str(tsi_no_na)


# Check number of samples per group (6 categories)
number_per_group <- tsi_no_na %>% 
  group_by(tsi_levels) %>%
  summarise(n = n())

# Counter plot in 7 categories
counter_plot <- number_per_group %>%
  ggplot(aes(x = factor(tsi_levels, levels = c("0-6", "6-12", 
                                               "12-24", "24-36",
                                               "36-48", "48-60",
                                               ">60"
                                              )
                       ), 
              y = n, 
              fill = factor(tsi_levels)
            )
        ) +
  geom_col(aes(fill = factor(tsi_levels))) +
  scale_fill_manual(label = c("0-6", "6-12", 
                              "12-24", "24-36",
                              "36-48", "48-60",
                              ">60"), 
                    name = "Subtypes", 
                    values = c("#0099CC", "#0D4F8B",
                               "#33A1C9", "#33A1DE", 
                               "#FFCC11", "#236B8E", 
                               "#33A1DE")) +
  labs(y = "Number of samples", x = "Months") +
  # coord_flip() +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 80),
                     breaks = seq(0,80,20)
                     ) +
  theme(text = element_text(size = 14, color = "grey15"),
        panel.background = element_rect(fill = "white",
                                        colour = "grey80"), 
        axis.title = element_text(size = 14, 
                                  color = "grey15"), 
        legend.text = element_text(size = 14),
        legend.position = "none") +
  geom_text(aes(label = n), vjust = -1, colour = "grey15", size = 5)

# Call a graph variable
counter_plot


# Boxplot in 7 categories
tsi_no_na %>%
  group_by(tsi_levels) %>%
  ggplot(aes(x = factor(tsi_levels, 
                        levels = c("0-6","6-12","12-24","24-36","36-48","48-60",">60")), 
             y = tsi_diff)) +
  geom_hline(yintercept = 0, size = 0.5, color = "#3A5894") +
  stat_boxplot(lwd = 1, geom = 'errorbar', 
               aes(color=factor(tsi_levels, 
                                levels = c("0-6","6-12","12-24","24-36","36-48","48-60",">60")))) +
  geom_boxplot(lwd=1, fatten=1, 
               aes(color=factor(tsi_levels, 
                                levels = c("0-6","6-12","12-24","24-36","36-48","48-60",">60"))), 
                                outlier.shape = 21) +
  scale_color_manual(label=c("0-6","6-12","12-24","24-36","36-48","48-60",">60"), 
                     name ="TSI Time", 
                     values = c("#0099CC", "#0D4F8B", "#33A1C9", "#33A1DE", 
                                "#FFCC11", "#236B8E", "#33A1DE")) +
  #geom_jitter(size=0.2, alpha=0.8, aes(color=factor(tsi_levels, levels = c("0-6","6-12","12-24","24-36","36-48","48-60", ">60")))) +
  geom_jitter(size = 0.15, alpha= 0.8, color = "grey15") +
  labs(x = 'Months', y = 'Known TSI - Estimated TSI\n(months)')+
  scale_y_continuous(limits = c(-90, 90), breaks = seq(-90,90,30)) +
  theme_minimal() +
  theme(text = element_text(size = 14, color = "grey15"), 
        panel.background = element_rect(fill = "white", 
                                        colour = "grey80"), 
        axis.title = element_text(size=14, color = "grey15"), 
        legend.text = element_text(size=14)) +
  theme(legend.position = "below")


# Contingency table: known vs estimated (threshold is 6 months)
tsi_contingency_6 <- tsi_df %>%
  rename(known_tsi = duration_of_infection_months) %>%
  rename(est_tsi = hivphylotsi_months) %>%
  select(known_tsi, est_tsi) %>%
  drop_na() %>%
  mutate(known_tsi_levels=ifelse(known_tsi<=6,"recent", "long-term")) %>%
  mutate(est_tsi_levels=ifelse(est_tsi<=6,"recent", "long-term"))


tsi_contingency_6$known_tsi_levels = paste("Known =", tsi_contingency_6$known_tsi_levels)
tsi_contingency_6$est_tsi_levels = paste("Estimated =", tsi_contingency_6$est_tsi_levels)
table(tsi_contingency_6$known_tsi_levels, tsi_contingency_6$est_tsi_levels)

# Accuracy (%) - 85.1%
# Accuracy = (true positives + true negatives) / (total)
tsi_accuracy_6 = (33 + 218)/ (33+218+8+36)*100
tsi_accuracy_6


# Contingency table: known vs estimated (threshold is 12 months)
tsi_contingency_12 <- tsi_df %>%
  rename(known_tsi = duration_of_infection_months) %>%
  rename(est_tsi = hivphylotsi_months) %>%
  select(known_tsi, est_tsi) %>%
  drop_na() %>%
  mutate(known_tsi_levels=ifelse(known_tsi<=12,"recent", "long-term")) %>%
  mutate(est_tsi_levels=ifelse(est_tsi<=12,"recent", "long-term"))


tsi_contingency_12$known_tsi_levels = paste("Known =", tsi_contingency_12$known_tsi_levels)
tsi_contingency_12$est_tsi_levels = paste("Estimated =", tsi_contingency_12$est_tsi_levels)
table(tsi_contingency_12$known_tsi_levels, tsi_contingency_12$est_tsi_levels)



# Accuracy (%) - 78.3%
# Accuracy = (true positives + true negatives) / (total)
tsi_accuracy_12 = (55 + 176)/ (55+176+31+33)*100
tsi_accuracy_12


# Contingency table: known vs estimated (threshold is 18 months)
tsi_contingency_18 <- tsi_df %>%
  rename(known_tsi = duration_of_infection_months) %>%
  rename(est_tsi = hivphylotsi_months) %>%
  select(known_tsi, est_tsi) %>%
  drop_na() %>%
  mutate(known_tsi_levels=ifelse(known_tsi<=18,"recent", "long-term")) %>%
  mutate(est_tsi_levels=ifelse(est_tsi<=18,"recent", "long-term"))


tsi_contingency_18$known_tsi_levels = paste("Known =", tsi_contingency_18$known_tsi_levels)
tsi_contingency_18$est_tsi_levels = paste("Estimated =", tsi_contingency_18$est_tsi_levels)
table(tsi_contingency_18$known_tsi_levels, tsi_contingency_18$est_tsi_levels)

# Accuracy (%) - 82.4%
# Accuracy = (true positives + true negatives) / (total)
tsi_accuracy_18 = (94 + 149)/ (94+149+23+29)*100
tsi_accuracy_18


# Contingency table: known vs estimated (threshold is 24 months)
tsi_contingency_24 <- tsi_df %>%
  rename(known_tsi = duration_of_infection_months) %>%
  rename(est_tsi = hivphylotsi_months) %>%
  select(known_tsi, est_tsi) %>%
  drop_na() %>%
  mutate(known_tsi_levels=ifelse(known_tsi<=24,"recent", "long-term")) %>%
  mutate(est_tsi_levels=ifelse(est_tsi<=24,"recent", "long-term"))


tsi_contingency_24$known_tsi_levels = paste("Known =", tsi_contingency_24$known_tsi_levels)
tsi_contingency_24$est_tsi_levels = paste("Estimated =", tsi_contingency_24$est_tsi_levels)
table(tsi_contingency_24$known_tsi_levels, tsi_contingency_24$est_tsi_levels)

# Accuracy (%) - 79.3%
# Accuracy = (true positives + true negatives) / (total)
tsi_accuracy_24 = (111 + 123)/ (111+123+38+23)*100
tsi_accuracy_24



# Contingency table: serology vs estimated (threshold is 6 months)
serology_contingency <- tsi_df %>%
  rename(known_tsi = duration_of_infection_months) %>%
  rename(est_tsi = hivphylotsi_months) %>%
  select(scount, known_tsi, est_tsi, serology) %>%
  drop_na() %>%
  filter(serology!="without") %>%
  mutate(known_tsi_levels=ifelse(known_tsi<=6,"recent", "long-term")) %>%
  mutate(est_tsi_levels=ifelse(est_tsi<=6,"recent", "long-term"))

serology_contingency$known_tsi_levels =paste("Known =", serology_contingency$known_tsi_levels)
serology_contingency$est_tsi_levels =paste("Estimated =", serology_contingency$est_tsi_levels)
serology_contingency$serology = paste("Serology =", serology_contingency$serology)
table(serology_contingency$known_tsi_levels, serology_contingency$est_tsi_levels)
table(serology_contingency$known_tsi_levels, serology_contingency$serology)

# Accuracy (%) known vs estimated (no NA, no "without") - 86%
# Accuracy = (true positives + true negatives) / (total)
tsi_accuracy_filtered = (27 + 212)/ (31+8+212+27)*100
tsi_accuracy_filtered

# Accuracy (%) known vs serology (no NA, no "without") - 82%
# Accuracy = (true positives + true negatives) / (total)
serology_accuracy_filtered = (40 + 188)/ (40+188+32+18)*100
serology_accuracy_filtered


# Subtypes in the dataset
subtypes_table <- tsi_df %>%
  rename(known_tsi = duration_of_infection_months) %>%
  rename(est_tsi = hivphylotsi_months) %>%
  select(known_tsi, est_tsi, subtype_prrt) %>%
  drop_na() %>%
  group_by(subtype_prrt)

subtype_plot <- subtypes_table %>%
  ggplot(aes(x=fct_infreq(subtype_prrt), fill=fct_infreq(subtype_prrt))) +
  geom_bar(aes(fill=fct_infreq(subtype_prrt))) +
  scale_fill_manual(name ="Subtypes", 
            values = c("#0099CC","#0D4F8B", "#53868B","#FFCC11","#236B8E",
                       "#33A1DE", "#008B8B","#8E8E38", "#33A1DE", "#5D92B1", "#8E8E88",
                       "#6E7B8B", "#33A1C9", "#4A708B", "#87CEFF", "#008B8B", "#8E8E38" )) +
  coord_flip() +
  labs(x="Subtype", y="Count") +
  scale_y_continuous(limits = c(0, 254)) +
  theme_minimal() +
  theme(text=element_text(size=14, color="grey15"),
        panel.background = element_rect(fill = "white",
        colour = "grey80"), 
        axis.title=element_text(size=14, color="grey15"), 
        legend.text=element_text(size=14),
        legend.position = "none")+
  geom_text(aes(label = ..count..), stat = "count", 
        vjust = 0.5, hjust=-0.3,
        colour = "grey15", size = 3) + scale_y_break(c(24,224))
subtype_plot


# Group 0_6
# Lolly plots for comparison, long table (by group)
long_tsi_0_6 <- tsi_no_na %>%
  filter(tsi_levels == "0-6") %>%
  select(scount, known_tsi, est_tsi) %>%
  pivot_longer(!scount, names_to = "tsi", values_to = "tsi_values") %>%
  mutate(tsi_values = ifelse(tsi == "est_tsi",
                               tsi_values, -1*tsi_values))

breaks_values <- pretty(long_tsi_0_6$tsi_values)


lollypop_long_0_6 <- long_tsi_0_6 %>%
  ggplot() +
  geom_segment(aes(x=scount, xend=scount, y=0, yend=tsi_values, color=tsi), size=0.5, linetype="longdash" ) +
  geom_point(aes(x=scount, y=tsi_values, color=tsi), size=1.5, shape=19) +
  geom_hline(yintercept = 6, size=0.5, color="#0D4F8B") +
  geom_hline(yintercept = -6, size=0.5, color="#FFCC11") +
  scale_y_continuous(breaks = breaks_values, limits = c(-15, 90),
                     labels = abs(breaks_values)) +
  scale_color_manual(name = 'TSI (group: 0-6 months)',
                     labels = c( "estimated","known"), 
                     values =c("#0D4F8B","#FFCC11")) +
  labs(y="TSI (months)", x="Sample") +
  theme_minimal() +
  theme(legend.position=c(0.75, 0.900), legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.background =element_rect(color="white"))+
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white",colour = "grey80"), 
        axis.title=element_text(size=14, color="grey15"), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.4, size=6))
lollypop_long_0_6



# Group 6_12
# Lolly plots for comparison, long table (by group)
long_tsi_6_12 <- tsi_no_na %>%
  filter(tsi_levels == "6-12") %>%
  select(scount, known_tsi, est_tsi) %>%
  pivot_longer(!scount, names_to = "tsi", values_to = "tsi_values") %>%
  mutate(tsi_values = ifelse(tsi == "est_tsi",
                             tsi_values, -1*tsi_values))

breaks_values <- pretty(long_tsi_6_12$tsi_values)


lollypop_long_6_12 <- long_tsi_6_12 %>%
  ggplot() +
  geom_segment(aes(x=scount, xend=scount, y=0, yend=tsi_values, color=tsi), size=0.5, linetype="longdash" ) +
  geom_point(aes(x=scount, y=tsi_values, color=tsi), size=1.5, shape=19) +
  geom_hline(yintercept = 6, size=0.5, color="#0D4F8B") +
  geom_hline(yintercept = 12, size=0.5, color="#0D4F8B") +
  geom_hline(yintercept = -6, size=0.5, color="#FFCC11") +
  geom_hline(yintercept = -12, size=0.5, color="#FFCC11") +
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values)) +
  scale_color_manual(name = 'TSI (group: 6-12 months)',
                     labels = c( "estimated","known"), 
                     values =c("#0D4F8B","#FFCC11")) +
  labs(y="TSI (months)", x="Sample") +
  theme_minimal() +
  theme(legend.position=c(0.4, 0.900), legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.background =element_rect(color="white"))+
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white",colour = "grey80"), 
        axis.title=element_text(size=14, color="grey15"), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.4, size=6))
lollypop_long_6_12




# Group 12-24
# Lolly plots for comparison, long table (by group)
long_tsi_12_24 <- tsi_no_na %>%
  filter(tsi_levels == "12-24") %>%
  select(scount, known_tsi, est_tsi) %>%
  pivot_longer(!scount, names_to = "tsi", values_to = "tsi_values") %>%
  mutate(tsi_values = ifelse(tsi == "est_tsi",
                             tsi_values, -1*tsi_values))

breaks_values <- pretty(long_tsi_12_24$tsi_values)

lollypop_long_12_24 <- long_tsi_12_24 %>%
  ggplot() +
  geom_segment(aes(x=scount, xend=scount, y=0, yend=tsi_values, color=tsi), size=0.5, linetype="longdash" ) +
  geom_point(aes(x=scount, y=tsi_values, color=tsi), size=1.5, shape=19) +
  geom_hline(yintercept = 12, size=0.5, color="#0D4F8B") +
  geom_hline(yintercept = 24, size=0.5, color="#0D4F8B") +
  geom_hline(yintercept = -12, size=0.5, color="#FFCC11") +
  geom_hline(yintercept = -24, size=0.5, color="#FFCC11") +
  scale_y_continuous(breaks = breaks_values,labels = abs(breaks_values)) +
  scale_color_manual(name = 'TSI (group: 12-24 months)',
                     labels = c( "estimated","known"), 
                     values =c("#0D4F8B","#FFCC11")) +
  labs(y="TSI (months)", x="Sample") +
  theme_minimal() +
  theme(legend.position=c(0.68, 0.935), legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.background =element_rect(color="white"))+
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white",colour = "grey80"), 
        axis.title=element_text(size=14, color="grey15"), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.4, size=6))
lollypop_long_12_24



scatter_6_12 <- tsi_no_na %>%
  filter(tsi_levels == "6-12") %>%
  ggplot() +
  geom_point(aes(x=known_tsi, y=est_tsi, color="blue"), size=1, shape=19) +
  labs(y="Estimated TSI", x="Known TSI") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0,70,20)) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0,70,20)) +
  theme_minimal() +
  theme(legend.position=c(0.4, 0.900), legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.background =element_rect(color="white"))+
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white",colour = "grey80"), 
        axis.title=element_text(size=14, color="grey15"), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.4, size=6))
scatter_6_12



# Correlation of Viral Load and TSI Difference   --------------------------


tsi_diff_virus_load <- tsi_df %>%
  rename(known_tsi = duration_of_infection_months) %>%
  rename(est_tsi = hivphylotsi_months) %>%
  select(scount, known_tsi, est_tsi, viral_load) %>%
  drop_na() %>%
  mutate(tsi_diff=known_tsi - est_tsi) %>%
  mutate(log10_viral_load = log10(viral_load)) %>%
  ggplot(aes(x = abs(tsi_diff), y = log10_viral_load)) + 
  geom_point(color="#0D4F8B") +
  geom_smooth(method='lm', formula= y~x, se=F, color="#3A5894")+
  stat_cor(label.y = 10)+ #this means at 10th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.y = 9.5) + #this means at 9.5th unit regression line equation will be shown
  labs(y="log10 Viral load", x="Absolute value of difference between Known TSI and Estimated TSI (months)") +
  theme_minimal() +
  theme(text=element_text(size=14, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"), 
        axis.title=element_text(size=14, color="grey15"), 
        legend.text=element_text(size=14),
        legend.position = "none")
tsi_diff_virus_load

