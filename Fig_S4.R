# -----------------------------
# 0. Load required libraries
# -----------------------------
library(phyloseq)   # Microbiome data manipulation
library(dplyr)      # Data wrangling
library(tidyr)      # Pivoting / reshaping data
library(ggplot2)    # Plotting
library(grid)       # For unit() in ggplot theme

# -----------------------------
# 1. Aggregate features at the desired taxonomic level
# -----------------------------
ESMARES_AAPs_tax <- tax_glom(ESMARES_AAPs_02.2026, taxrank = "NewTax")

# -----------------------------
# 2. Transform counts to relative abundance (%) per sample
# -----------------------------
ESMARES_AAPs_perc <- transform_sample_counts(
  ESMARES_AAPs_tax,
  function(x) (x / sum(x)) * 100
)

# -----------------------------
# 3. Convert ASV table to tidy (long) format
# -----------------------------
otu_df <- as.data.frame(otu_table(ESMARES_AAPs_perc))
otu_df_t <- as.data.frame(t(otu_df))       # transpose so samples are rows
otu_df_t$ASV <- rownames(otu_df_t)         # add ASV column

tidy_df <- gather(
  otu_df_t,
  key = "Sample",
  value = "Abundance",
  -ASV
)
tidy_df <- tidy_df[, c("Sample", "ASV", "Abundance")]

# -----------------------------
# 4. Add taxonomic information
# -----------------------------
tax_df <- as.data.frame(tax_table(ESMARES_AAPs_perc))
tax_df$ASV <- rownames(tax_df)

tidy_df2 <- merge(
  tidy_df,
  tax_df[, c("ASV", "NewTax")],
  by = "ASV",
  all.x = TRUE
)

# -----------------------------
# 5. Add sample metadata
# -----------------------------
sample_df <- as.data.frame(sample_data(ESMARES_AAPs_perc))
sample_df$Sample <- rownames(sample_df)

sample_df2 <- sample_df[, c("Sample", "Community", "Station", "Cruise", "Sample2")]

barplot_df <- left_join(tidy_df2, sample_df2, by = "Sample")

# -----------------------------
# 6. Factor ordering for plotting
# -----------------------------
# Order cruises
barplot_df$Cruise <- factor(
  barplot_df$Cruise,
  levels = c("ES0120","ES0620","ES0121","ES0521","ES0122","ES0322","ES1122")
)

# Order stations
barplot_df$Station <- factor(
  barplot_df$Station,
  levels = c("AG1","AG2","AG3","AG4","ST1","ST2","ST3","ST5","ST6","ST7",
             "MA1","MA2","MA3","MA4","MA5","MA7","AL1","AL2","AL3","AL5")
)

# Order taxa
barplot_df$NewTax <- factor(
  barplot_df$NewTax,
  levels = c("A","B","C","D","G1 (Rhodobacterales)","G2 (Rhodobacterales)",
             "Other Rhodobacterales","I (Burkholderiales)","J",
             "K (Pseudomonadales)","Rhizobiales","Sphingomonadales",
             "Unclassified","Others")
)

# Define bar colors
bar_colors <- c("#0F8299FF", "#3E9FB3FF", "#7ABECCFF", "#B8DEE6FF", "#7B4173FF",
                "#A55194FF", "#CE6DBDFF", "#D6616BFF","#F47942FF","#FBB04EFF",
                "#B88100FF","#00A070FF","#E0E0D0FF","#707058FF")

# Order Sample2 factor (controls order on y-axis)
barplot_df$Sample2 <- factor(
  barplot_df$Sample2,
  levels = sort(unique(barplot_df$Sample2), decreasing = TRUE)
)

# -----------------------------
# 7. Combine metadata columns for y-axis labels
# -----------------------------
# For example: Cruise | Station
barplot_df$y_label <- paste(barplot_df$Cruise, barplot_df$Station, sep = " | ")

# Create named vector to preserve Sample2 order while changing displayed labels
y_labels <- setNames(barplot_df$y_label, barplot_df$Sample2)

# -----------------------------
# 8. Plot stacked barplot of relative abundance
# -----------------------------
ggplot(barplot_df, aes(x = Abundance, y = Sample2, color = NewTax, fill = NewTax)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Community, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 10) +
  ylab(NULL)+
  xlab("Relative Abundance (%)") +
  theme(
    axis.text.y = element_text(size = 5, angle = 0, vjust = 0.4),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    panel.spacing = unit(2, "cm")
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_color_manual(values = bar_colors) +
  guides(fill = guide_legend(ncol = 10)) +
  scale_y_discrete(labels = y_labels)  # Show combined metadata labels on y-axis

