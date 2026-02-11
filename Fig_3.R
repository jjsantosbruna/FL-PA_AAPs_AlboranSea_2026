# ============================================================
# LIBRARIES
# ============================================================
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
library(paletteer)
library(patchwork)

# ============================================================
# PREPROCESSING PHYLOSEQ OBJECT
# ============================================================

# Remove problematic sample
ESMARES_AAPs_filtered <- prune_samples(
  !(sample_names(ESMARES_AAPs_02.2026) %in% "Esmares_0121_208"),
  ESMARES_AAPs_02.2026
)

# Rarefaction (reproducible)
set.seed(23)
ESMARES_AAPs_rarefied <- rarefy_even_depth(ESMARES_AAPs_filtered)

# ============================================================
# SPLIT COMMUNITIES
# ============================================================

FL <- subset_samples(ESMARES_AAPs_rarefied, Community == "Free-living")
PA <- subset_samples(ESMARES_AAPs_rarefied, Community == "Particle-attached")

# ============================================================
# NMDS — FREE-LIVING COMMUNITY
# ============================================================

# ---------- Prepare OTU table ----------
otu_FL <- data.frame(otu_table(FL))
otu_FL$Code <- rownames(otu_FL)

# ---------- Prepare metadata ----------
metadata_FL <- data.frame(sample_data(FL))
metadata_FL$Code <- rownames(metadata_FL)

# ---------- Merge OTUs and metadata ----------
beta_table_FL <- merge(metadata_FL, otu_FL, by = "Code")

# ---------- Separate environmental variables and abundances ----------
env_FL   <- beta_table_FL[, 12:37]
abund_FL <- beta_table_FL[, 38:5274]

# Remove unwanted environmental variables
env_FL <- env_FL %>%
  select(-APPs_conc, -AAPs_per)

# ---------- NMDS ----------
# Square-root normalization
abund_FL_sqrt <- sqrt(abund_FL)

# Bray–Curtis dissimilarity
dist_FL <- vegdist(abund_FL_sqrt, method = "bray", na.rm = TRUE)

# NMDS ordination
nmds_FL <- metaMDS(dist_FL)

# ---------- Environmental fitting ----------
env_FL_scaled <- as.data.frame(scale(env_FL))
envfit_FL <- envfit(nmds_FL, env_FL_scaled, na.rm = TRUE)

# ---------- Extract NMDS coordinates ----------
nmds_scores_FL <- as.data.frame(nmds_FL$points)
beta_table_FL$nmds_x <- nmds_scores_FL$MDS1
beta_table_FL$nmds_y <- nmds_scores_FL$MDS2

# ---------- Extract significant environmental vectors ----------
env_vectors_FL <- as.data.frame(scores(envfit_FL, display = "vectors"))
env_vectors_FL$p_value <- envfit_FL$vectors$pvals

env_vectors_sig_FL <- env_vectors_FL %>%
  filter(p_value < 0.05) %>%
  mutate(variable = rownames(.))

# ---------- Aesthetics ----------
color_palette_FL <- c(
  "#172554", "#2563EB", "#60A5FA",
  "#A5D8FF", "#006D2CFF", "#74C476FF", "#D68910"
)

cruise_shapes_FL <- c(16, 17, 15, 18, 16, 17, 16, 17)

# ---------- Plot NMDS FL ----------
A=ggplot(beta_table_FL, aes(nmds_x, nmds_y, color = Cruise)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(shape = Cruise)) +
  geom_text(aes(label = Station),
            vjust = -0.8, size = 2, show.legend = FALSE) +
  geom_segment(
    data = env_vectors_sig_FL,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "gray50"
  ) +
  geom_text(
    data = env_vectors_sig_FL,
    aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = variable),
    color = "gray50", size = 3
  ) +
  theme_bw() +
  scale_color_manual(values = color_palette_FL) +
  scale_shape_manual(values = cruise_shapes_FL) +
  theme(legend.position = "right") +
  labs(x = "NMDS1", y = "NMDS2", title="Free-living")

# ============================================================
# NMDS — PARTICLE-ATTACHED COMMUNITY
# ============================================================

# ---------- Prepare OTU table ----------
otu_PA <- data.frame(otu_table(PA))
otu_PA$Code <- rownames(otu_PA)

# ---------- Prepare metadata ----------
metadata_PA <- data.frame(sample_data(PA))
metadata_PA$Code <- rownames(metadata_PA)

# ---------- Merge OTUs and metadata ----------
beta_table_PA <- merge(metadata_PA, otu_PA, by = "Code")

# ---------- Separate environmental variables and abundances ----------
env_PA   <- beta_table_PA[, 12:37]
abund_PA <- beta_table_PA[, 38:5274]

# Remove unwanted variables
env_PA <- env_PA %>%
  select(-AAPs_per, -APPs_conc)

# ---------- NMDS ----------
abund_PA_sqrt <- sqrt(abund_PA)
dist_PA <- vegdist(abund_PA_sqrt, method = "bray", na.rm = TRUE)
nmds_PA <- metaMDS(dist_PA)

# ---------- Environmental fitting ----------
env_PA_scaled <- as.data.frame(scale(env_PA))
envfit_PA <- envfit(nmds_PA, env_PA_scaled, na.rm = TRUE)

# ---------- Extract NMDS coordinates ----------
nmds_scores_PA <- as.data.frame(nmds_PA$points)
beta_table_PA$nmds_x <- nmds_scores_PA$MDS1
beta_table_PA$nmds_y <- nmds_scores_PA$MDS2

# ---------- Extract significant environmental vectors ----------
env_vectors_PA <- as.data.frame(scores(envfit_PA, display = "vectors"))
env_vectors_PA$p_value <- envfit_PA$vectors$pvals

env_vectors_sig_PA <- env_vectors_PA %>%
  filter(p_value < 0.05) %>%
  mutate(variable = rownames(.))

# ---------- Aesthetics ----------
color_palette_PA <- c(
  "#2563EB", "#60A5FA", "#A5D8FF",
  "#006D2CFF", "#D68910"
)

cruise_shapes_PA <- c(17, 15, 18, 16, 16)

# ---------- Plot NMDS PA ----------
B=ggplot(beta_table_PA, aes(nmds_x, nmds_y, color = Cruise)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(shape = Cruise)) +
  geom_text(aes(label = Station),
            vjust = -0.8, size = 2, show.legend = FALSE) +
  geom_segment(
    data = env_vectors_sig_PA,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "gray50"
  ) +
  geom_text(
    data = env_vectors_sig_PA,
    aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = variable),
    color = "gray50", size = 3
  ) +
  theme_bw() +
  scale_color_manual(values = color_palette_PA) +
  scale_shape_manual(values = cruise_shapes_PA) +
  theme(legend.position = "right") +
  labs(x = "NMDS1", y = "NMDS2", title = "Particle-associated")

# ---------- Arrange both plots in a row ----------

  A + B
