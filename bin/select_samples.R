#!/usr/bin/env Rscript

# 1. Libraries ------------------------------------------------------------

suppressMessages({
  library(argparse)
  library(readxl)
  library(dplyr)
})



# 2. Create parser --------------------------------------------------------

parser <- ArgumentParser(description = "Filter Excel data to select samples 
                                        for TSI retraining output CSV + TXT")

parser$add_argument("-i", "--input", required = TRUE, 
                                     help = "Input Excel file (.xlsx)")
parser$add_argument("-s", "--sheet", default = "1", 
                                     help = "Excel sheet name or number (default: 1)")
parser$add_argument("-c", "--output-csv", required = TRUE, 
                                          help = "Output filtered CSV file")
parser$add_argument("-t", "--output-txt", required = TRUE, 
                                          help = "Output TXT file with scount values")


args <- parser$parse_args()

# Convert sheet to numeric if possible
sheet_arg <- suppressWarnings(as.numeric(args$sheet))
if (is.na(sheet_arg)) sheet_arg <- args$sheet


# 3. Define the filtering function ----------------------------------------

filter_excel_data <- function(input_file, output_csv, output_txt, sheet = 1) {
  df <- read_excel(input_file, sheet = sheet)
 
  cols_to_select <- c("scount", "project", "known_tsi_days",
                      "known_tsi_months", "est_tsi_months", 
                      "known_tsi_years", "sample_number", 
                      "first_or_followup","viral_load", 
                      "subtype_prrt", "subtype_kallisto",
                      "multiplicity", "multi_all_genes", 
                      "elisa_sum_bed", "biorad_avidity",
                      "protocol", "file_name",
                      "file_size_mb_r1", "file_size_mb_r2",
                      "reads_raw_mln", "reads_final_mln",  
                      "pol_gap", "pol_gap_info", "pol_cov", 
                      "gag_gap", "gag_gap_info", "gag_cov",
                      "comment")
  
  df_selected <- df |>  select(all_of(cols_to_select))
  
  
  # Exclude samples with low coverage or with gaps in either pol or gag regions
  # coverage by phyloscanner ( > 31.62 ) 
  # multiplicity by phyloscanner ( visual inspections, selected if not all genes involved )
  
  df_filtered <- df_selected |> 
    filter(multi_all_genes != "yes") |> 
    filter(pol_gap != "yes", pol_cov != "uneven", gag_gap != "yes", gag_cov != "uneven") |> 
    filter(known_tsi_days != 0, known_tsi_years != 1) |> # noisy samples
    arrange(scount)
  
  write.csv(df_filtered, output_csv, row.names = FALSE)
  writeLines(df_filtered$scount, con = output_txt)
}


# 4. Call the function with parsed arguments ------------------------------
filter_excel_data(args$input, args$output_csv, args$output_txt, sheet_arg)


