
# Load required libraries for plotting and plot composition
library(ggplot2)
library(patchwork)

# ----------------------------------
# Panel A: Violin plot of DNA concentration (QUBIT_DNA)
# grouped by amplification success (Amplified_pufM)
# ----------------------------------

A = ggplot(Table_S4, aes(x = Amplified_pufM, y = QUBIT_DNA)) +
  theme_bw() +  # Use a clean black-and-white theme
  geom_violin(scale = "width") +  # Violin plot scaled to the same width
  geom_jitter(
    aes(color = Filter_Size_(µm)), # Color points by filter size
    size = 2,
    shape = 1
  ) +
  scale_color_manual(
    values = c("0.2" = "#F47942FF", "3" = "#00A070FF") # Custom colors per filter size
  ) +
  scale_y_continuous(
    trans = "asinh",  # Apply transformation to visualize better from 1 to 10 ng/uL
    breaks = c(0, 1, 2, 4, 8, 16, 32, 64) # Define y-axis breaks
  )

# ----------------------------------
# Subset data to only include samples
# filtered with 3 µm filters (PA)
# ----------------------------------

data_PA <- subset(DNA_pufM_table, Filter_Size == "3")

# ----------------------------------
# Panel B: Violin plot of AAPs abundance (cells per mL)
# for the 3 µm filter subset
# ----------------------------------

B = ggplot(data_PA, aes(x = Amplified_pufM, y = AAPs_cells_mL)) +
  theme_bw() +
  geom_violin(scale = "width") +
  geom_jitter(
    aes(color = Filter_Size),  # Color points by filter size
    size = 2,
    shape = 1
  ) +
  scale_color_manual(
    values = c("0.2" = "#F47942FF", "3" = "#00A070FF")
  )

# ----------------------------------
# Panel C: Violin plot of AAPs percentage
# for the 3 µm filter subset
# ----------------------------------

C = ggplot(data_PA, aes(x = Amplified_pufM, y = AAPs_per)) +
  theme_bw() +
  geom_violin(scale = "width") +
  geom_jitter(
    aes(color = Filter_Size),
    size = 2,
    shape = 1
  ) +
  scale_color_manual(
    values = c("0.2" = "#F47942FF", "3" = "#00A070FF")
  ) +
  scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15) # Custom y-axis breaks
  )

# ----------------------------------
# Combine panels A, B, and C into a single row
# ----------------------------------

A | B | C

# ----------------------------------
# Non-parametric statistical tests
# (Kruskal-Wallis tests)
# ----------------------------------

# Test differences in DNA concentration across amplification success (all data)
kruskal.test(QUBIT_DNA ~ Amplified_pufM, data = DNA_pufM_table)

# Test differences in DNA concentration across amplification success (PA only)
kruskal.test(QUBIT_DNA ~ Amplified_pufM, data = data_PA)

# Test differences in AAPs concentration (cells/mL) across Amplified_pufM (PA only)
kruskal.test(AAPs_con ~ Amplified_pufM, data = data_PA)

# Test differences in AAPs percentage across mplification success (PA only)
kruskal.test(AAPs_per ~ Amplified_pufM, data = data_PA)
