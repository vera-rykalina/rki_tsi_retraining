# Load library
library(TidyMultiqc)
library(readxl)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(tidyr)
library(ggsci)

# Export original table
tsi_df <- read_excel("data/hiv_phylotsi_validation.xlsx",sheet = "HIVtime_v2")


df = TidyMultiqc::load_multiqc("data/multiqc_data.json")

df <- df |>
  separate_wider_delim(metadata.sample_id, delim = "_", names = c("id", "operation"))

df <- spread(df, key = operation, value = general.total_sequences)

