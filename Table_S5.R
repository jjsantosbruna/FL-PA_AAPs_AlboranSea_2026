library(phyloseq)
library(dplyr)
library(tidyr)

# Define the taxonomic level to be used
tax <- "NewTax"

# 1. Agglomerate features at the selected taxonomic level
ps_tax <- tax_glom(ESMARES_AAPs_02.2026, taxrank = tax)

# 2. Transform counts to relative abundance PER SAMPLE
ps_rel <- transform_sample_counts(ps_tax, function(x) x / sum(x))

# 3. Convert to data frame
df <- psmelt(ps_rel)
# Abundance = relative abundance of each ASV in each sample

# Keep only FL and PA (optional but recommended)
df <- df %>%
  filter(Community %in% c("Free-living", "Particle-attached"))

# 4. SUM relative abundances per taxon WITHIN each sample
df_tax_sample <- df %>%
  group_by(Sample, Community, !!sym(tax)) %>%
  summarise(
    Abundance_tax_sample = sum(Abundance),
    .groups = "drop"
  )
# Now each row = one taxon in one sample

# 5. Mean relative abundance across ALL samples
global <- df_tax_sample %>%
  group_by(!!sym(tax)) %>%
  summarise(Global = mean(Abundance_tax_sample), .groups = "drop") %>%
  mutate(Global = Global * 100)

# 6. Mean relative abundance by fraction (FL / PA)
by_fraction <- df_tax_sample %>%
  group_by(Community, !!sym(tax)) %>%
  summarise(Mean_abundance = mean(Abundance_tax_sample), .groups = "drop") %>%
  mutate(Mean_abundance = Mean_abundance * 100) %>%
  pivot_wider(
    names_from = Community,
    values_from = Mean_abundance,
    values_fill = 0
  )

# 7. Kruskal–Wallis test (FL vs PA) using per-sample taxon abundances
pvalues <- df_tax_sample %>%
  group_by(!!sym(tax)) %>%
  summarise(
    p_value = tryCatch(
      kruskal.test(Abundance_tax_sample ~ Community)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# 8. Merge everything into the final table
table_S4 <- global %>%
  left_join(by_fraction, by = tax) %>%
  left_join(pvalues, by = tax)
