# GenomOrthoANICalculator

A lightweight Python-based command-line pipeline implementing the OrthoANI algorithm for calculating pairwise Average Nucleotide Identity (ANI) among multiple genome sequences.

This implementation is based on:

Lee et al. (2016). OrthoANI: An improved algorithm and software for calculating average nucleotide identity. 
International Journal of Systematic and Evolutionary Microbiology.DOI: 10.1099/ijsem.0.000760

In contrast to the original Java-based OAT software described by Lee et al. (2016), this implementation provides automated batch processing of multiple genome FASTA files within a single directory.

---

## Features

- 1020 bp genome fragmentation
- Reciprocal BLASTn best-hit strategy
- 35% minimum alignment length cutoff (357 bp)
- Directional ANI (A→B and B→A)
- Reciprocal OrthoANI calculation
- Supports unlimited genome comparisons
- Command-line based (Unix/Mac/Linux compatible)

---

## Requirements

- Python ≥ 3.7
- pandas
- Biopython
- NCBI BLAST+ (makeblastdb, blastn available in PATH)

Install Python dependencies:

pip install pandas biopython

---

## Input

A directory containing genome FASTA files:

Supported formats:
- .fasta
- .fna
- .fa

---

## Usage

python GenomOrthoANICalculator.py \
    --input /path/to/genomes \
    --output /path/to/output \
    --threads 4

---

## Output

Results file:

OrthoANI_results.tsv

Columns:
- Genome1
- Genome2
- ANI_1_to_2
- ANI_2_to_1
- Reciprocal_ANI

Each genome pair also has its own working directory containing:
- Fragmented sequences
- BLAST databases
- BLAST output files

---

## Algorithm Overview

1. Each genome is fragmented into 1020 bp consecutive segments.
2. Fragments shorter than 1020 bp are discarded.
3. BLASTn is performed using:
   - task = blastn
   - dust = no
   - xdrop_gap = 150
   - penalty = -1
   - reward = 1
   - evalue = 1e-15
4. Alignments ≥ 35% of fragment length (≥357 bp) are retained.
5. Reciprocal best hits (RBH) are identified.
6. ANI is calculated as the mean nucleotide identity of orthologous fragment pairs.

---

## Visualization Scripts

Additional scripts are provided in the `scripts/` directory to generate ANI matrices and heatmaps.

#- create_ani_matrix.py → Converts OrthoANI_results.tsv into square ANI matrix
- heatmap_python.py → Generates clustered heatmap using seaborn/matplotlib
- heatmap_R.R → Generates heatmap using pheatmap in R
- Diagonal_Heatmap_dendrogram.py
   (python Diagonal_Heatmap_dendrogram.py OrthoANI_results.tsv)
These scripts are optional and intended for visualization of genome similarity.


## Citation

If you use this tool, please cite:

Lee et al. (2016) OrthoANI: An improved algorithm and software for calculating average nucleotide identity.
Int J Syst Evol Microbiol. DOI: 10.1099/ijsem.0.000760

And cite this software as:

Umar, M. GenomOrthoANICalculator v1.0 (2026). 
Available at: https://github.com/MdUmar-tech/GenomOrthoANICalculator

---

## License

MIT License

## Background

This tool was originally developed during my PhD after encountering several issue due to java while using the Java-based GUI OAT (OrthoANI Tool) on macOS-2020 (Intel).

To enable reproducible, large-scale, command-line ANI analysis, I independently implemented the OrthoANI algorithm (Lee et al., 2016) in Python with full batch automation support.

Unlike the original GUI-based OAT and web-based implementations (e.g., EZBioCloud OrthoANI), this version supports unlimited genome comparisons and is designed for Unix/Linux server and HPC environments.
