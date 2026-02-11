# ===============================================================
# Load required libraries
# ===============================================================
library(phyloseq)    
library(tidyverse)  
library(pheatmap)   
library(dplyr)      
library(Hmisc)       
library(ggpubr)      
library(gridExtra)   
library(grid)        

# ===============================================================
# 1. Prepare phyloseq object with desired taxonomy
# ===============================================================

# Extract the taxonomy table from the phyloseq object
tax_mat <- as.matrix(tax_table(ESMARES_AAPs_02.2026))

# Keep only the column of interest ("NewTax")
tax_mat <- tax_mat[, "NewTax", drop = FALSE]

# Create a temporary dummy column to allow tax_glom() to work
tax_mat <- cbind(tax_mat, Dummy = "dummy")

# Reassign the updated taxonomy table back to the phyloseq object
tax_table(ESMARES_AAPs_02.2026) <- tax_table(tax_mat)

# Agglomerate OTUs at the "NewTax" level
phy_tax <- tax_glom(ESMARES_AAPs_02.2026, taxrank = "NewTax")

# ===============================================================
# 2. Process Free-living community
# ===============================================================

# Subset samples labeled as "Free-living"
FL <- subset_samples(phy_tax, Community == "Free-living")

# Agglomerate at the "NewTax" level
ESMARES_AAPs_mod_tax_FL <- tax_glom(FL, taxrank = "NewTax")

# Extract OTU table and convert to a data frame
otu_family_FL <- otu_table(ESMARES_AAPs_mod_tax_FL)
df_family_FL <- as.data.frame(otu_family_FL)

# Rename columns with taxonomy names
colnames(df_family_FL) <- as.data.frame(tax_table(ESMARES_AAPs_mod_tax_FL))$NewTax

# Remove low abundance taxa
df_family_FL <- df_family_FL[, !names(df_family_FL) %in% c(
  "Unclassified","Sphingomonadales","I (Burkholderiales)","J","Others","Rhizobiales")]

# Compute relative abundances as percentages
df_per_FL <- df_family_FL / rowSums(df_family_FL) * 100

# ===============================================================
# 3. Prepare metadata
# ===============================================================
metadata_FL <- as.data.frame(sample_data(FL))
metadata2_FL <- metadata_FL[, 11:36]  # Select relevant columns
metadata2_FL$Coast_dis <- metadata_FL$Coast_nm  # Add coastline distance

# ===============================================================
# 4. Compute Spearman correlations
# ===============================================================
correlations_FL <- rcorr(as.matrix(df_per_FL), as.matrix(metadata2_FL), type = "spearman")

# Extract correlation coefficients and p-values
r_matrix_FL <- correlations_FL[["r"]][9:35, 1:8]
p_matrix_FL <- correlations_FL[["P"]][9:35, 1:8]

# Apply Bonferroni correction for multiple testing
p_flat_FL <- as.vector(p_matrix_FL)
p_bon_FL <- p.adjust(p_flat_FL, method = "bonferroni")
p_matrix_FL2 <- matrix(p_bon_FL, nrow = nrow(p_matrix_FL))

# ===============================================================
# 5. Generate significance asterisks
# ===============================================================
get_asterisks_FL <- function(p) {
  if (is.na(p)) return("")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

# Create a matrix of asterisks for significance
asterisks_matrix_FL <- matrix("", nrow = nrow(p_matrix_FL), ncol = ncol(p_matrix_FL))
for (i in 1:nrow(p_matrix_FL)) {
  for (j in 1:ncol(p_matrix_FL)) {
    asterisks_matrix_FL[i, j] <- get_asterisks_FL(p_matrix_FL[i, j])
  }
}

# Transpose matrices for plotting
asterisks_matrix_t_FL <- t(asterisks_matrix_FL)
r_matrix_t_FL <- t(r_matrix_FL)

# Remove not significant variables from the matrix
r_matrix_t_FL <- r_matrix_t_FL[, !colnames(r_matrix_t_FL) %in% c("E4E6","E2E3","S275_295")]
asterisks_matrix_t_FL <- asterisks_matrix_t_FL[, -c(3,4,5)]

# ===============================================================
# 6. Plot heatmap
# ===============================================================
A <- grid.grabExpr(
  pheatmap(
    r_matrix_t_FL,
    border_color = "black",
    number_color = "white",
    cluster_rows = TRUE,           # Cluster metadata variables
    cluster_cols = FALSE,          # Keep taxa order fixed
    display_numbers = asterisks_matrix_t_FL,  # Show significance
    fontsize_number = 10,
    color = colorRampPalette(c("#1F0CA9", "white", "#E82D6E"))(50),  # Blue-white-pink gradient
    labels_row = colnames(df_per_FL),
    labels_col = colnames(r_matrix_t_FL),
    legend = FALSE,
    main = "Free-living"
  )
)


# ===============================================================
# PARTICLE-ATTACHED COMMUNITY
# ===============================================================

# Subset Particle-associated samples 
PA <- subset_samples(phy_tax, Community == "Particle-attached")

# Extract taxonomy table for reference
tax_table_met_PA <- as.data.frame(tax_table(PA))

# Agglomerate OTUs at the "NewTax" level
ESMARES_AAPs_mod_tax_PA <- tax_glom(PA, taxrank = "NewTax")

# Extract OTU table and convert to a data frame
otu_family_PA <- otu_table(ESMARES_AAPs_mod_tax_PA)
df_family_PA <- as.data.frame(otu_family_PA)

# Rename columns with taxonomy names
colnames(df_family_PA) <- as.data.frame(tax_table(ESMARES_AAPs_mod_tax_PA))$NewTax

# Remove low abundance taxa
df_family_PA <- df_family_PA[, !names(df_family_PA) %in% c(
  "Unclassified","Sphingomonadales","I (Burkholderiales)","J","Others","Rhizobiales")]

# Compute relative abundances as percentages
df_per_PA <- df_family_PA / rowSums(df_family_PA) * 100

# ===============================================================
# Prepare metadata
# ===============================================================

metadata_PA <- as.data.frame(sample_data(PA))

# Remove AAP columns due to not existing data
metadata_PA <- metadata_PA[, !names(metadata_PA) %in% c("AAPs_per","APPs_conc")]

# Select columns for correlation analysis
metadata2_PA <- metadata_PA[, 11:34]

# Add coastline distance variable
metadata2_PA$Coast_dis <- metadata_PA$Coast_nm

# ===============================================================
# Compute Spearman correlations
# ===============================================================
correlations_PA <- rcorr(as.matrix(df_per_PA), as.matrix(metadata2_PA), type = "spearman")

# Extract correlation coefficients and p-values
r_matrix_PA <- correlations_PA[["r"]][9:33, 1:8]
p_matrix_PA <- correlations_PA[["P"]][9:33, 1:8]

# ===============================================================
# Apply Bonferroni correction to p-values
# ===============================================================

# Flatten matrix, adjust p-values, and reconstruct matrix
p_flat_PA <- as.vector(p_matrix_PA)
p_bon_PA <- p.adjust(p_flat_PA, method = "bonferroni")
p_matrix_PA <- matrix(p_bon_PA, nrow = nrow(p_matrix_PA))

# Set column and row names (optional)
colnames(p_matrix_PA) <- rownames(p_matrix_PA) <- colnames(p_flat_PA)

# ===============================================================
# Generate significance asterisks
# ===============================================================
get_asterisks_PA <- function(p) {
  if (is.na(p)) return("")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

# Apply asterisks function to entire p-value matrix
asterisks_matrix_PA <- matrix("", nrow = nrow(p_matrix_PA), ncol = ncol(p_matrix_PA))
for (i in 1:nrow(p_matrix_PA)) {
  for (j in 1:ncol(p_matrix_PA)) {
    asterisks_matrix_PA[i, j] <- get_asterisks_PA(p_matrix_PA[i, j])
  }
}

# Transpose matrices for plotting
asterisks_matrix_t_PA <- t(asterisks_matrix_PA)
r_matrix_t_PA <- t(r_matrix_PA)

# Remove not significant variables from the matrix
r_matrix_t_PA <- r_matrix_t_PA[, !colnames(r_matrix_t_PA) %in% c(
  "Picoeuk","Temp","PO4","NO3","SR","NO2","Chla_T","Bact_T",
  "Chla_20u","a325","a254","E4E6","E2E3","S275_295","SiO3")]

asterisks_matrix_t_PA <- asterisks_matrix_t_PA[, -c(1,2,3,4,5,7,8,9,12,13,14,15,17,18,22)]

# ===============================================================
# Plot heatmap for Particle-attached community
# ===============================================================
B <- grid.grabExpr(
  pheatmap(
    r_matrix_t_PA,
    border_color = "black",
    number_color = "white",
    cluster_rows = TRUE,           # Cluster metadata variables
    cluster_cols = FALSE,          # Keep taxa order fixed
    display_numbers = asterisks_matrix_t_PA, # Show significance
    fontsize_number = 10,
    color = colorRampPalette(c("#1F0CA9", "white", "#E82D6E"))(50), # Blue-white-pink gradient
    labels_row = colnames(df_per_PA),
    labels_col = colnames(asterisks_matrix_t_PA),
    legend_breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
    main = "Particle-associated"
  )
)

# ===============================================================
# Combine Free-living and Particle-attached heatmaps vertically
# ===============================================================
grid.arrange(A, B, ncol = 1, heights = c(1, 1))

