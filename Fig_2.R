# ============================================================
# 0. Load required libraries
# ============================================================
library(phyloseq)   # Microbiome data manipulation
library(ggplot2)    # Data visualization
library(vegan)      # Diversity indices and Bray–Curtis distances
library(car)        # Additional statistical utilities (optional)
library(patchwork)  # Combine multiple ggplots into one figure
library(FSA)        # Dunn post-hoc test for non-parametric comparisons


# ============================================================
# PART 1: ALPHA DIVERSITY ANALYSIS
# ============================================================

# ------------------------------------------------------------
# 1. Remove problematic or low-quality samples (if applicable)
# ------------------------------------------------------------
ESMARES_AAPs_filt <- prune_samples(
  !(sample_names(ESMARES_AAPs_02.2026) %in% "Esmares_0121_208"),
  ESMARES_AAPs_02.2026
)

# ------------------------------------------------------------
# 2. Rarefy samples to an even sequencing depth
#    (minimum sequencing depth across samples)
# ------------------------------------------------------------
ESMARES_AAPs_rare <- rarefy_even_depth(
  ESMARES_AAPs_filt,
  sample.size = min(sample_sums(ESMARES_AAPs_filt)),
  rngseed = 23,        # Seed for reproducibility
  replace = TRUE,
  trimOTUs = TRUE,
  verbose = TRUE
)

# ------------------------------------------------------------
# 3. Extract OTU table and sample metadata
# ------------------------------------------------------------
otu_mat <- as.data.frame(otu_table(ESMARES_AAPs_rare))

sample_df <- as.data.frame(sample_data(ESMARES_AAPs_rare))
sample_df$SampleID <- rownames(sample_df)

# ------------------------------------------------------------
# 4. Calculate alpha diversity indices per sample
#    - Shannon diversity index
#    - Chao1 richness estimator
# ------------------------------------------------------------
alpha_df <- data.frame(
  SampleID = rownames(otu_mat),
  Shannon  = diversity(otu_mat, index = "shannon"),
  Chao1    = specnumber(otu_mat) +
    (rowSums(otu_mat == 1) * (rowSums(otu_mat == 1) - 1)) /
    (2 * (specnumber(otu_mat) + 1))
)

# Merge diversity indices with sample metadata
alpha <- merge(alpha_df, sample_df, by = "SampleID")

# ------------------------------------------------------------
# 5. Plot Chao1 richness by community (FL vs PA)
# ------------------------------------------------------------
A <- ggplot(alpha, aes(x = Community, y = Chao1, fill = Community)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_summary(fun = "mean", geom = "point", color = "white", size = 2) +
  scale_fill_manual(values = c("#F47942FF", "#00A070FF")) +
  labs(x = NULL, y = "Chao1 index") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Free-living" = "FL",
    "Particle-attached" = "PA"
  ))

# ------------------------------------------------------------
# 6. Plot Shannon diversity by community (FL vs PA)
# ------------------------------------------------------------
B <- ggplot(alpha, aes(x = Community, y = Shannon, fill = Community)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_summary(fun = "mean", geom = "point", color = "white", size = 2) +
  scale_fill_manual(values = c("#F47942FF", "#00A070FF")) +
  labs(x = NULL, y = "Shannon index") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Free-living" = "FL",
    "Particle-attached" = "PA"
  ))


# ============================================================
# PART 2: BETA DIVERSITY ANALYSIS (BRAY–CURTIS)
# ============================================================

# ------------------------------------------------------------
# 7. Subset samples by community
# ------------------------------------------------------------
FL <- subset_samples(ESMARES_AAPs_02.2026, Community == "Free-living")
PA <- subset_samples(ESMARES_AAPs_02.2026, Community == "Particle-attached")

# Extract OTU tables for each subset
otu_FL  <- as.data.frame(otu_table(FL))
otu_PA  <- as.data.frame(otu_table(PA))
otu_all <- as.data.frame(otu_table(ESMARES_AAPs_02.2026))

# ------------------------------------------------------------
# 8. Compute Bray–Curtis dissimilarities
#    - Within FL
#    - Within PA
#    - Between FL and PA
# ------------------------------------------------------------
bray_FL <- vegdist(otu_FL, method = "bray")
bray_PA <- vegdist(otu_PA, method = "bray")

# Compute Bray–Curtis across all samples
bray_all <- vegdist(otu_all, method = "bray")
bray_mat <- as.matrix(bray_all)

# Extract between-community distances (FL vs PA)
bray_FL_PA <- as.vector(bray_mat[rownames(otu_FL), rownames(otu_PA)])

# ------------------------------------------------------------
# 9. Combine all distance values into a single data frame
# ------------------------------------------------------------
bray_df <- data.frame(
  Distance = c(bray_FL, bray_PA, bray_FL_PA),
  Group = factor(
    rep(c("FL", "PA", "FL&PA"),
        times = c(length(bray_FL),
                  length(bray_PA),
                  length(bray_FL_PA)))
  )
)

# ------------------------------------------------------------
# 10. Plot Bray–Curtis distance distributions
# ------------------------------------------------------------
C <- ggplot(bray_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_summary(fun = "mean", geom = "point", color = "white", size = 2) +
  scale_fill_manual(values = c("#F47942FF", "#B070F4", "#00A070FF")) +
  labs(x = NULL, y = "Bray–Curtis distance") +
  theme_bw() +
  theme(legend.position = "none")


# ============================================================
# PART 3: COMBINE FIGURES
# ============================================================

# Combine alpha and beta diversity plots into a single row
A + B + C


# ============================================================
# PART 4: STATISTICAL TESTS
# ============================================================

# ------------------------------------------------------------
# 11. Kruskal–Wallis tests
#    (non-parametric comparison between groups)
# ------------------------------------------------------------
kruskal.test(Chao1 ~ Community, data = alpha)
kruskal.test(Shannon ~ Community, data = alpha)

kruskal.test(Distance ~ Group, data = bray_df)

# ------------------------------------------------------------
# 12. Dunn post-hoc test (pairwise comparisons)
#    Applied after a significant Kruskal–Wallis test
# ------------------------------------------------------------
dunnTest(Distance ~ Group, data = bray_df)

