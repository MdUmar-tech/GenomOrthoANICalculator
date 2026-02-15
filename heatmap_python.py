import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load file
df = pd.read_csv("results/OrthoANI_results.tsv", sep="\t")

# Pivot using Reciprocal_ANI
ani_matrix = df.pivot(index="Genome1",
                      columns="Genome2",
                      values="Reciprocal_ANI")

# Make symmetric
ani_matrix = ani_matrix.combine_first(ani_matrix.T)

# Ensure matrix is float and writable
ani_matrix = ani_matrix.astype(float).copy()

# Set diagonal = 100
for genome in ani_matrix.index:
    ani_matrix.loc[genome, genome] = 100

# Sort rows & columns same order
ani_matrix = ani_matrix.sort_index(axis=0)
ani_matrix = ani_matrix.sort_index(axis=1)

# Save matrix
ani_matrix.to_csv("results/ANI_matrix.tsv", sep="\t")

print("Matrix created successfully.")

# Create heatmap
plt.figure(figsize=(20, 15))
sns.heatmap(
    ani_matrix,
    cmap="viridis",
    annot=True,
    fmt=".1f",
    linewidths=0.5
)

plt.title("ANI Heatmap (Reciprocal ANI)", fontsize=18)
plt.xticks(rotation=90)
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig("results/ani_heatmap.png", dpi=300)
plt.show()

print("Heatmap saved as results/ani_heatmap.png")
