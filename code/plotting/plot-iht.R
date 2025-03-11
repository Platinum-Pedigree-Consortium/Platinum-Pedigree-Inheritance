#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(viridis)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if enough arguments are provided
if (length(args) < 2) {
  stop("Usage: Rscript plot_script.R <input_file> <output_file>")
}

file_path <- args[1]      # First argument: input iht file 
output_file <- args[2]    # Second argument: output PDF file

# Load the data
data <- read.table(file_path, header=TRUE, comment.char="")

# Remove the last column
data <- data[, -ncol(data)]

# Convert start and end to numeric
data <- data %>% mutate(across(c(start, end), as.numeric))

# Reshape the data for plotting, splitting haplotypes
long_data <- data %>% 
  pivot_longer(cols = -c(X.chrom, start, end, marker_count, len), names_to = "Individual", values_to = "Genotype") %>% 
  separate(Genotype, into = c("Haplotype1", "Haplotype2"), sep = "[|/]") %>% 
  pivot_longer(cols = c(Haplotype1, Haplotype2), names_to = "Haplotype", values_to = "Allele") %>% 
  mutate(Haplotype = factor(Haplotype, levels = c("Haplotype1", "Haplotype2")))

# Convert Allele to a categorical factor
allele_levels <- unique(long_data$Allele)
long_data$Allele <- factor(long_data$Allele, levels = allele_levels)

# Order individuals so that both haplotypes are grouped together
long_data <- long_data %>% 
  arrange(Individual, Haplotype) %>% 
  mutate(Individual_Hap = factor(paste(Individual, Haplotype, sep = "_"), levels = unique(paste(Individual, Haplotype, sep = "_"))))

# Create the plot
p <- ggplot(long_data, aes(y = Individual_Hap, color = Allele, group = Individual_Hap)) +
  geom_segment(aes(x = start/1e6, xend = end/1e6, y = Individual_Hap, yend = Individual_Hap), size = 2) +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "", x = "Chrom. pos. MB", y = "Individual (Haplotypes)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 8))

# Save the plot as an image
ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
message("Plot saved as ", output_file)
