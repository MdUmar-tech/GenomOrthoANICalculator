# Install packages if needed
# install.packages("tidyverse")
# install.packages("pheatmap")

library(tidyverse)
library(pheatmap)

# Load file
df <- read.delim("results/OrthoANI_results.tsv", header = TRUE)

# Create pivot matrix using Reciprocal_ANI
ani_matrix <- df %>%
  select(Genome1, Genome2, Reciprocal_ANI) %>%
  pivot_wider(names_from = Genome2,
              values_from = Reciprocal_ANI)

# Convert to matrix
ani_matrix <- as.data.frame(ani_matrix)
rownames(ani_matrix) <- ani_matrix$Genome1
ani_matrix$Genome1 <- NULL

ani_matrix <- as.matrix(ani_matrix)

# Make symmetric
ani_matrix[is.na(ani_matrix)] <- t(ani_matrix)[is.na(ani_matrix)]

# Ensure numeric
ani_matrix <- apply(ani_matrix, 2, as.numeric)
rownames(ani_matrix) <- colnames(ani_matrix)

# Set diagonal = 100
diag(ani_matrix) <- 100

# Sort rows and columns
ani_matrix <- ani_matrix[order(rownames(ani_matrix)),
                         order(colnames(ani_matrix))]

# Save matrix
write.table(ani_matrix,
            file = "results/ANI_matrix.tsv",
            sep = "\t",
            quote = FALSE)

cat("Matrix created successfully\n")

# Heatmap
pheatmap(
  ani_matrix,
  #color = colorRampPalette(c("blue", "yellow", "red"))(100),
  display_numbers = TRUE,
  number_format = "%.1f",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  filename = "results/ani_heatmap.png",
  width = 20,
  height = 15
)

cat("Heatmap saved as results/ani_heatmap.png\n")
