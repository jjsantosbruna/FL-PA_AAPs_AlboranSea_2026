# Load required libraries
library(phyloseq)
library(ggplot2)
library(ggpubr)

# Extract metadata from phyloseq object
metadata <- as(sample_data(ESMARES_AAPs_02.2026), "data.frame")

# Set cruise order as factor
metadata$Cruise <- factor(metadata$Cruise,
                          levels = c("ES0120","ES0620","ES0121",
                                     "ES0521","ES0122","ES0322","ES1122"))

# -------------------------
# Boxplots for environmental variables
# -------------------------

# Sea Surface Temperature
A <- ggplot(metadata, aes(x = Cruise, y = Temp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL, y = "SST (ºC)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Sea Surface Salinity
B <- ggplot(metadata, aes(x = Cruise, y = Sal)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL, y = "SSS") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Total Chlorophyll-a
C <- ggplot(metadata, aes(x = Cruise, y = Chla_T)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("Total Chl ", italic(a), " (", mu, "g ", L^-1, ")"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Dissolved Oxygen
D <- ggplot(metadata, aes(x = Cruise, y = O2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("DO (", mu, "M )"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Nitrate
E <- ggplot(metadata, aes(x = Cruise, y = NO3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("NO"[3], " (", mu, "M)"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Nitrite
F <- ggplot(metadata, aes(x = Cruise, y = NO2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("NO"[2], " (", mu, "M)"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Ammonium
G <- ggplot(metadata, aes(x = Cruise, y = NH4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("NH"[4], " (", mu, "M)"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Phosphate (corrected variable assumed to be PO4)
H <- ggplot(metadata, aes(x = Cruise, y = PO4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("PO"[4], " (", mu, "M)"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Silicate
I <- ggplot(metadata, aes(x = Cruise, y = SiO3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3,
               color = "#019875FF", position = position_dodge(width = 0.75)) +
  labs(x = NULL,
       y = expression(paste("SiOH (", mu, "M)"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange all plots in a 3x3 grid
ggarrange(A, B, C, D, E, F, G, H, I,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
          ncol = 3, nrow = 3)