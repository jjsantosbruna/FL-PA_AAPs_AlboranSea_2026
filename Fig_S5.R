# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

##########################
# 1. Extract OTU table
##########################

otu <- as.data.frame(otu_table(ESMARES_AAPs_02.2026))

# Transpose if taxa are rows
if (taxa_are_rows(ESMARES_AAPs_02.2026)) {
  otu <- t(otu)
}

otu <- as.data.frame(otu)
otu$SampleID <- rownames(otu)  # Add SampleID column

##########################
# 2. Extract metadata
##########################

metadata <- as.data.frame(sample_data(ESMARES_AAPs_02.2026))
metadata$SampleID <- rownames(metadata)

##########################
# 3. Convert to long format
##########################

otu_long <- otu %>%
  pivot_longer(-SampleID, names_to = "ASV", values_to = "Counts") %>%
  left_join(metadata, by = "SampleID")

##########################
# 4. Bubble plot by Cruise
##########################

# Total abundance per ASV per Cruise
asv_by_cruise <- otu_long %>%
  group_by(Cruise, ASV) %>%
  dplyr::summarise(Abundance = sum(Counts), .groups = "drop")

# Total abundance per Cruise
total_by_cruise <- asv_by_cruise %>%
  group_by(Cruise) %>%
  dplyr::summarise(Total = sum(Abundance), .groups = "drop")

# Calculate relative abundance (%)
asv_relative <- asv_by_cruise %>%
  left_join(total_by_cruise, by = "Cruise") %>%
  mutate(RelativeAbundance = (Abundance / Total) * 100)

# Select top 10 ASVs globally
top10_asvs <- asv_relative %>%
  group_by(ASV) %>%
  dplyr::summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:10) %>%
  pull(ASV)

# Filter top 10
asv_relative_top10 <- asv_relative %>%
  filter(ASV %in% top10_asvs)

# Add taxonomy
taxonomy <- as.data.frame(tax_table(ESMARES_AAPs_02.2026))
taxonomy$ASV <- rownames(taxonomy)

asv_relative_top10 <- asv_relative_top10 %>%
  left_join(taxonomy, by = "ASV")

# Order ASVs
asv_order <- c("asv10", "asv8", "asv7", "asv6", "asv4",
               "asv14", "asv5", "asv3", "asv2", "asv1")
asv_relative_top10$ASV <- factor(asv_relative_top10$ASV, levels = asv_order)

# Order Cruise
cruise_order <- c("ES0120", "ES0620", "ES0121", "ES0521", 
                  "ES0122", "ES0322", "ES1122")
asv_relative_top10$Cruise <- factor(asv_relative_top10$Cruise, levels = cruise_order)

# Colors for taxonomy (adjust as needed)
bar_colors <- c("#0F8299", "#3E9FB3", "#7ABECC", "#7B4173", "#FBB04E")

# Bubble plot for Cruise
plot_cruise <- ggplot(asv_relative_top10, 
                      aes(x = Cruise, y = ASV, size = RelativeAbundance, color = NewTax)) +
  geom_point(alpha = 0.8) +
  scale_size_area(max_size = 10, limits = c(0, 25)) +
  scale_color_manual(values = bar_colors) +
  labs(x = NULL, y = NULL,
       size = "Relative Abundance (%)",
       color = "Taxonomy",
       title = "Cruise") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")

##########################
# 5. Bubble plot by Transect
##########################

# Total abundance per ASV per Transect
asv_by_transect <- otu_long %>%
  group_by(Transect, ASV) %>%
  dplyr::summarise(Abundance = sum(Counts), .groups = "drop")

# Total abundance per Transect
total_by_transect <- asv_by_transect %>%
  group_by(Transect) %>%
  dplyr::summarise(Total = sum(Abundance), .groups = "drop")

# Relative abundance (%)
asv_relative <- asv_by_transect %>%
  left_join(total_by_transect, by = "Transect") %>%
  mutate(RelativeAbundance = (Abundance / Total) * 100)

# Top 10 ASVs
top10_asvs <- asv_relative %>%
  group_by(ASV) %>%
  dplyr::summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:10) %>%
  pull(ASV)

asv_relative_top10 <- asv_relative %>%
  filter(ASV %in% top10_asvs) %>%
  left_join(taxonomy, by = "ASV")

# Order ASVs
asv_relative_top10$ASV <- factor(asv_relative_top10$ASV, levels = asv_order)

# Order Transect
transect_order <- c("Algeciras", "Sotogrande", "Malaga", "Almeria")
asv_relative_top10$Transect <- factor(asv_relative_top10$Transect, levels = transect_order)

# Bubble plot for Transect
plot_transect <- ggplot(asv_relative_top10, 
                        aes(x = Transect, y = ASV, size = RelativeAbundance, color = NewTax)) +
  geom_point(alpha = 0.8) +
  scale_size_area(max_size = 10, limits = c(0, 25)) +
  scale_color_manual(values = bar_colors) +
  labs(x = NULL, y = NULL,
       size = "Relative Abundance (%)",
       color = "Taxonomy",
       title = "Transect") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

##########################
# 6. Bubble plot by Community
##########################

# Total abundance per ASV per Community
asv_by_comm <- otu_long %>%
  group_by(Community, ASV) %>%
  dplyr::summarise(Abundance = sum(Counts), .groups = "drop")

# Total abundance per Community
total_by_comm <- asv_by_comm %>%
  group_by(Community) %>%
  dplyr::summarise(Total = sum(Abundance), .groups = "drop")

# Relative abundance (%)
asv_relative <- asv_by_comm %>%
  left_join(total_by_comm, by = "Community") %>%
  mutate(RelativeAbundance = (Abundance / Total) * 100)

# Top 10 ASVs
top10_asvs <- asv_relative %>%
  group_by(ASV) %>%
  dplyr::summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:10) %>%
  pull(ASV)

asv_relative_top10 <- asv_relative %>%
  filter(ASV %in% top10_asvs) %>%
  left_join(taxonomy, by = "ASV")

# Order ASVs
asv_relative_top10$ASV <- factor(asv_relative_top10$ASV, levels = asv_order)

# Bubble plot for Community
plot_community <- ggplot(asv_relative_top10, 
                         aes(x = Community, y = ASV, size = RelativeAbundance, color = NewTax)) +
  geom_point(alpha = 0.8) +
  scale_size_area(max_size = 10, limits = c(0, 25)) +
  scale_color_manual(values = bar_colors) +
  labs(x = NULL, y = NULL,
       size = "Relative Abundance (%)",
       color = "Taxonomy",
       title = "Community") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

##########################
# 7. Combine plots with patchwork
##########################

plot_cruise + plot_transect + plot_community + 
  plot_layout(widths = c(2, 1, 0.5))  # Adjust relative widths