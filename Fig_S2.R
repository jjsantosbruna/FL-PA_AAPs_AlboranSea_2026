# Load required libraries
library(phyloseq)
library(ggplot2)
library(ggpubr)

# Extract metadata from phyloseq object
metadata <- as(sample_data(ESMARES_AAPs_02.2026), "data.frame")

# Set cruise order as factor to control x-axis order
metadata$Cruise <- factor(metadata$Cruise,
                          levels = c("ES0120","ES0620","ES0121",
                                     "ES0521","ES0122","ES0322","ES1122"))

# -------------------------
# Boxplots for microbial abundance variables
# -------------------------

# AAP abundance (cells mL⁻¹)
A <- ggplot(metadata, aes(x = Cruise, y = APPs_conc)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("AAPs (cells ", mL^-1, ")"))) +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# AAP relative abundance (%)
B <- ggplot(metadata, aes(x = Cruise, y = AAPs_per)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = "AAPs (%)") +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Total bacterial abundance (cells mL⁻¹)
C <- ggplot(metadata, aes(x = Cruise, y = Bact_T)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("Total Bacteria (cells ", mL^-1, ")"))) +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Low Nucleic Acid bacteria (LNA)
D <- ggplot(metadata, aes(x = Cruise, y = Bact_LNA)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("Bacteria LNA (cells ", mL^-1, ")"))) +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# High Nucleic Acid bacteria (HNA)
E <- ggplot(metadata, aes(x = Cruise, y = Bact_HNA)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("Bacteria HNA (cells ", mL^-1, ")"))) +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Prochlorococcus abundance
F <- ggplot(metadata, aes(x = Cruise, y = Prochlo)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste(italic("Prochlorococcus"),
                            " (cells ", mL^-1, ")"))) +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Synechococcus abundance
G <- ggplot(metadata, aes(x = Cruise, y = Synecho)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste(italic("Synechococcus"),
                            " (cells ", mL^-1, ")"))) +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Picoeukaryotes abundance (PPEs)
H <- ggplot(metadata, aes(x = Cruise, y = Picoeuk)) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF",
               position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("Picoeukaryotes (cells ", mL^-1, ")"))) +
  theme_bw() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange all plots in a 2x4 grid
ggarrange(A, B, C, D, E, F, G, H,
          labels = c("A", "B", "C", "D",
                     "E", "F", "G", "H"),
          ncol = 2, nrow = 4)