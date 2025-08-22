

# 1. Libraries ------------------------------------------------------------
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


#font_import(prompt = FALSE) 
#loadfonts()

# 2. Export dataframe -----------------------------------------------------
hivtype_df <- read_excel("data/MS115_report.xlsx")
head(hivtype_df) 
nrow(hivtype_df)


subtype_plot <- hivtype_df |> 
  ggplot(aes(x = fct_infreq(Subtype), fill = fct_infreq(Subtype))) +
  geom_bar(aes(fill = fct_infreq(Subtype))) +
  coord_flip() +
  labs(x = "Subtype", y = "Count") +
  scale_fill_manual(name = "Subtypes", 
                   values = c("#0099CC","#0D4F8B", "#53868B","#FFCC11","#236B8E",
                              "#33A1DE", "#008B8B","#8E8E38", "#33A1DE", "#5D92B1", "#8E8E88",
                              "#6E7B8B", "#33A1C9", "#4A708B", "#87CEFF", "#008B8B", "#8E8E38" )) +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0,60,15)) +
  theme_minimal() +
  ggtitle("HIV-1 Subtyping (Stanford, Comet, Rega, Geno2Pheno)") +
  theme(text = element_text(size  = 12, color = "grey15"),
        panel.background = element_rect(fill  = "white", color = "grey80"), 
        axis.title  = element_text(size = 12, color = "grey15"), 
        legend.text = element_text(size = 12),
        legend.position = "none") +
  geom_text(aes(label = ..count..), 
            stat = "count", 
            vjust = 0.5, hjust = -0.3,
            colour = "grey15", size = 3) + 
  scale_y_break(c(15,55))
subtype_plot

ggsave("images/img_png/hivtype_countplot_poster_700_dpi.png", 
       units = "cm", width = 36, height = 18, dpi = 700)
