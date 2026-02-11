library(phyloseq)
library(EcolUtils)
library(ggplot2)
library(paletteer)
library(patchwork)

########################################
# 1. Extract OTU table and metadata
########################################

# Extract OTU table and sample metadata from phyloseq object
otu_table_df <- as.data.frame(otu_table(ESMARES_AAPs_02.2026))
metadata_df  <- as.data.frame(sample_data(ESMARES_AAPs_02.2026))

########################################
# 2. PAN index (Particle-Attached Niche)
########################################

# Extract Community column as a vector
community_vector <- metadata_df$Community

# Convert Community to numeric vector (1 = Particle-attached, 0 = Free-living)
pan_binary_vector <- as.numeric(community_vector == "Particle-attached")

# Convert OTU table to matrix (required by niche.val)
otu_matrix <- as.matrix(otu_table_df)

# Calculate niche values using permutations
set.seed(23)
pan_index_table <- niche.val(
  otu_matrix,
  pan_binary_vector,
  n = 1000,
  probs = c(0.05, 0.95)
)

# Keep only significant ASVs
pan_index_filtered <- subset(pan_index_table, sign != "NON SIGNIFICANT")
pan_index_filtered$ASV <- rownames(pan_index_filtered)

########################################
# 3. Merge PAN index with taxonomy
########################################

# Extract taxonomy table
taxonomy_df <- as.data.frame(tax_table(ESMARES_AAPs_02.2026))
taxonomy_df$ASV <- rownames(taxonomy_df)

# Merge taxonomy with PAN index results
tax_pan <- merge(taxonomy_df, pan_index_filtered, by = "ASV", all.y = TRUE)

# Remove uncultured MAGs
tax_pan <- subset(tax_pan, NewTax != "Uncultured MAGs")

# Define taxonomic order for plotting
tax_pan$NewTax <- factor(
  tax_pan$NewTax,
  levels = c(
    "Others", "Unclassified", "Sphingomonadales", "Rhizobiales",
    "K (Pseudomonadales)", "J", "I (Burkholderiales)",
    "Other Rhodobacterales", "G2 (Rhodobacterales)",
    "G1 (Rhodobacterales)", "D", "C", "B", "A"
  )
)

# Define color palette for taxonomic groups
tax_colors <- c(
  "#0F8299FF", "#3E9FB3FF", "#7ABECCFF", "#B8DEE6FF",
  "#7B4173FF", "#A55194FF", "#CE6DBDFF", "#D6616BFF",
  "#F47942FF", "#FBB04EFF", "#B88100FF", "#00A070FF",
  "#E0E0D0FF", "#707058FF"
)

# Violin plot for PAN index
plot_pan <- ggplot(
  tax_pan,
  aes(x = observed, y = NewTax, fill = forcats::fct_rev(NewTax))
) +
  geom_violin(scale = "width") +
  stat_summary(fun = median, geom = "point", size = 2, color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
  scale_fill_manual(values = tax_colors) +
  labs(x = "PAN index", y = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

########################################
# 4. CON index (Coast Niche)
########################################

# Extract coast index vector
coast_vector <- metadata_df$Coast_nm

# Calculate niche values for coast index
set.seed(23)
coast_index_table <- niche.val(
  otu_matrix,
  coast_vector,
  n = 1000,
  probs = c(0.05, 0.95)
)

# Keep only significant ASVs
coast_index_filtered <- subset(coast_index_table, sign != "NON SIGNIFICANT")
coast_index_filtered$ASV <- rownames(coast_index_filtered)

########################################
# 5. Merge CON index with taxonomy
########################################

tax_coast <- merge(taxonomy_df, coast_index_filtered, by = "ASV", all.y = TRUE)
tax_coast <- subset(tax_coast, NewTax != "Uncultured MAGs")

tax_coast$NewTax <- factor(
  tax_coast$NewTax,
  levels = c(
    "Others", "Unclassified", "Sphingomonadales", "Rhizobiales",
    "K (Pseudomonadales)", "J", "I (Burkholderiales)",
    "Other Rhodobacterales", "G2 (Rhodobacterales)",
    "G1 (Rhodobacterales)", "D", "C", "B", "A"
  )
)

# Violin plot for CON index
plot_con <- ggplot(
  tax_coast,
  aes(x = observed, y = NewTax, fill = forcats::fct_rev(NewTax))
) +
  geom_violin(scale = "width") +
  stat_summary(fun = median, geom = "point", size = 1, color = "white") +
  scale_fill_manual(values = tax_colors) +
  scale_x_reverse() +
  labs(x = "CON index", y = "", fill = "Taxonomy") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

########################################
# 6. Combine plots
########################################

(plot_pan + plot_con) +
  plot_annotation(tag_levels = "A")




