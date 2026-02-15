# Install packages if needed
# install.packages("tidyverse")
# install.packages("pheatmap")

library(tidyverse)
library(pheatmap)

# Load file
df <- read.delim("OrthoANI_results.tsv", header = TRUE)

df$Genome1 <- trimws(df$Genome1)
df$Genome2 <- trimws(df$Genome2)

# Keep required columns
df2 <- df %>%
  select(Genome1, Genome2, Reciprocal_ANI)

# Add reverse combinations
df_rev <- df2 %>%
  rename(Genome1 = Genome2,
         Genome2 = Genome1)

df_full <- bind_rows(df2, df_rev)

# Get all genomes
genomes <- unique(c(df_full$Genome1, df_full$Genome2))

# Add diagonal
diag_df <- data.frame(
  Genome1 = genomes,
  Genome2 = genomes,
  Reciprocal_ANI = 100
)

df_full <- bind_rows(df_full, diag_df)

# Remove duplicates
df_full <- df_full %>%
  distinct(Genome1, Genome2, .keep_all = TRUE)

# Create matrix
ani_matrix <- df_full %>%
  pivot_wider(names_from = Genome2,
              values_from = Reciprocal_ANI)

ani_matrix <- as.data.frame(ani_matrix)
rownames(ani_matrix) <- ani_matrix$Genome1
ani_matrix$Genome1 <- NULL

ani_matrix <- as.matrix(ani_matrix)

# Replace remaining NA safely
ani_matrix[is.na(ani_matrix)] <- 100

# Sort properly (optional)
ani_matrix <- ani_matrix[order(rownames(ani_matrix)),
                         order(colnames(ani_matrix))]

# Save matrix
write.table(ani_matrix,
            file = "ANI_matrix.tsv",
            sep = "\t",
            quote = FALSE)

cat("Matrix created successfully\n")

# Heatmap
pheatmap(
  ani_matrix,
  display_numbers = TRUE,
  number_format = "%.1f",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  filename = "ani_heatmap.png",
  width = 13,
  height = 10
)

cat("Heatmap saved successfully\n")
