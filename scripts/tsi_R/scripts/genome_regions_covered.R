library(tidyverse)
library(IRanges)
library(scales)

# Load annotation
ref <- read_csv("data/HXB2_refdata.csv")

genome_length <- max(ref$`HXB2 base position`)

# Colors for genes
custom_colors <- c(
  "gag" = "#0D4F8B",
  "pol" = "#33A1C9",
  "gp120" = "#FFCC11",
  "gp41" = "#33A1DE"
)

# Prepare gene positions from ref
gene_df <- ref %>%
  select(pos = `HXB2 base position`, RF1 = `RF1 protein`, RF2 = `RF2 protein`, RF3 = `RF3 protein`) %>%
  pivot_longer(cols = starts_with("RF"), names_to = "frame", values_to = "gene") %>%
  mutate(pos = as.integer(pos), gene = tolower(gene)) %>%
  filter(gene %in% names(custom_colors)) %>%
  distinct(pos, gene)

# Create continuous gene regions
gene_regions <- gene_df %>%
  arrange(gene, pos) %>%
  mutate(
    gap = (pos != lag(pos) + 1) | (gene != lag(gene)),
    gap = if_else(is.na(gap), TRUE, gap),
    group = cumsum(gap)
  ) %>%
  group_by(group, gene) %>%
  summarize(start = min(pos), end = max(pos), .groups = "drop") %>%
  mutate(label = gene)

# Separate gag and pol for overlapping visualization
gag_df <- gene_regions %>% filter(gene == "gag")
pol_df <- gene_regions %>% filter(gene == "pol")
other_genes_df <- gene_regions %>% filter(!gene %in% c("gag", "pol"))

combined_start <- min(gag_df$start, pol_df$start)
combined_end <- max(gag_df$end, pol_df$end)

overlap_start <- max(gag_df$start, pol_df$start)
overlap_end <- min(gag_df$end, pol_df$end)
overlap_df <- tibble(start = overlap_start, end = overlap_end) %>% filter(end >= start)

blend_colours <- function(cols, weights = NULL) {
  rgb_mat <- t(col2rgb(cols)) / 255
  if (is.null(weights)) weights <- rep(1 / nrow(rgb_mat), nrow(rgb_mat))
  mixed <- colSums(rgb_mat * weights)
  rgb(mixed[1], mixed[2], mixed[3])
}
gag_col <- custom_colors["gag"]
pol_col <- custom_colors["pol"]
overlap_col <- alpha(blend_colours(c(gag_col, pol_col)), 0.5)

# Load primer & product regions
regions <- read_csv("data/regions_igs_no_gag.csv") %>%
  separate(sense_primer_coords, into = c("sp_start", "sp_end"), sep = "-", convert = TRUE) %>%
  separate(antisense_primer_coords, into = c("asp_start", "asp_end"), sep = "-", convert = TRUE) %>%
  separate(product_coords_excl_primers, into = c("prod_start", "prod_end"), sep = "-", convert = TRUE) %>%
  mutate(group = "IGS samples, partial genomes")

# IRanges for primers
primer_ir <- IRanges(
  start = c(regions$sp_start, regions$asp_start),
  end = c(regions$sp_end, regions$asp_end)
)

# Reduce primers to merge any overlapping primer ranges
primer_ir_reduced <- reduce(primer_ir)

# IRanges for products (excluding primers in the CSV)
product_ir <- IRanges(start = regions$prod_start, end = regions$prod_end)

# Subtract primer positions from product ranges (get product without primers = coverage with gaps)
product_no_primers_ir <- setdiff(product_ir, primer_ir_reduced)

# Gene IRanges for clipping coverage strictly inside gene boundaries
gene_ir <- IRanges(start = gene_regions$start, end = gene_regions$end, names = gene_regions$gene)

# Now clip each coverage subregion (product_no_primers_ir) to corresponding gene(s)
clipped_list <- list()
for (i in seq_along(product_no_primers_ir)) {
  overlaps <- findOverlaps(product_no_primers_ir[i], gene_ir)
  if (length(overlaps) == 0) next
  for (j in seq_along(overlaps)) {
    gene_idx <- subjectHits(overlaps)[j]
    clipped_start <- max(start(product_no_primers_ir[i]), start(gene_ir)[gene_idx])
    clipped_end <- min(end(product_no_primers_ir[i]), end(gene_ir)[gene_idx])
    if (clipped_end >= clipped_start) {
      clipped_list[[length(clipped_list) + 1]] <- tibble(
        start = clipped_start,
        end = clipped_end,
        gene = names(gene_ir)[gene_idx],
        group = "IGS samples, partial genomes"
      )
    }
  }
}
clipped_coverage_df <- bind_rows(clipped_list)

# Calculate coverage percentage (exclude primers, clipped inside genes)
total_covered <- sum(clipped_coverage_df$end - clipped_coverage_df$start + 1)
coverage_pct <- (total_covered / genome_length) * 100
cat(sprintf("Coverage excluding primers and clipped to gene boundaries: %.2f%%\n", coverage_pct))

# Plot
p <- ggplot() +
  geom_segment(aes(x = 1, xend = genome_length, y = -0.15, yend = -0.15),
               color = "gray80", size = 3) +
  
  # Gene boxes for gag and pol without individual outlines, only combined outline
  geom_rect(data = gag_df, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.075),
            fill = gag_col, color = NA, size = 0) +  # No outline for gag
  geom_rect(data = pol_df, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.075),
            fill = pol_col, color = NA, size = 0) +  # No outline for pol
  
  # Blended color region for overlap between gag and pol, no outline here
  geom_rect(aes(xmin = overlap_start, xmax = overlap_end, ymin = 0, ymax = 0.075),
            fill = overlap_col, color = NA, size = 0) +  # No outline for overlap region
  
  # One continuous frame for the combined gag-pol gene box (single outline)
  geom_rect(aes(xmin = combined_start, xmax = combined_end, ymin = 0, ymax = 0.075),
            fill = NA, color = "black", size = 0.4) +  # Single outline for the whole gag-pol box
  
  # Other genes (non-gag, non-pol)
  geom_rect(data = other_genes_df, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.075, fill = gene),
            color = "black", size = 0.3, alpha = 0.8) +
  
  # Gene names
  geom_text(data = gene_regions, aes(x = (start + end) / 2, y = 0.0375, label = label),
            color = "black", size = 5, fontface = "italic", hjust = 0.5, vjust = 0.5) +
  
  # Region boxes with no shading (solid color, no alpha)
  geom_rect(data = clipped_coverage_df, aes(xmin = start, xmax = end, ymin = -0.1, ymax = -0.06),
            fill = "#778899", color = "black", size = 0.3, alpha = 1) +  # Solid color, no alpha blending
  
  # Coverage label
  geom_label(aes(x = genome_length * 0.5, y = 0.25,
                 label = sprintf("Coverage (IGS primer settings): %.2f%% vs 69.99%%", coverage_pct)),
             label.padding = unit(0.55, "lines"), label.size = 0.35,
             color = "grey10", fill = "white", size = 5, hjust = 0.5) +
  
  # Axes and plot labels
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05)),
                     limits = c(0, genome_length),
                     breaks = c(0, 2500, 5000, 7500, 10000)) +
  
  scale_fill_manual(values = custom_colors, name = NULL) +
  
  ylim(-0.17, 0.35) +
  
  labs(title = "HIV-1 Genome (HXB2, K03455): Genes and Primer-Excluded Coverage Regions",
       x = "HIV-1 genome position", y = NULL) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 16, color = "grey10"),
    panel.background = element_rect(fill = "white", colour = "grey20"),
    axis.title = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_line(color = "grey10"),
    axis.text.x = element_text(color = "grey10")
  )

# Print the plot
print(p)

# Save the plot
ggsave("images/img_png/hiv_genes_coverage_clipped_igs_no_gag.png", plot = p,
       width = 14, height = 3.4, dpi = 500)

