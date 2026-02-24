#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO

__version__ = "1.0"
# --------------------------------------------------
# Argument Parser
# --------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Pairwise OrthoANI calculation for all genomes in a folder"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input folder containing .fasta/.fna files"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output folder"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=2,
        help="Number of BLAST threads (default=2)"
    )
    
    parser.add_argument(
    "--version",
    action="version",
    version="GenomOrthoANICalculator v1.0.0"
   )
    return parser.parse_args()


# --------------------------------------------------
# Utility Functions
# --------------------------------------------------

def get_genome_files(folder):
    extensions = ["*.fasta", "*.fna", "*.fa"]
    files = []

    for ext in extensions:
        files.extend(glob.glob(os.path.join(folder, ext)))

    return sorted(files)


def create_structure(main_dir):
    sub_dirs = ["fragment", "mkdb", "blast"]
    os.makedirs(main_dir, exist_ok=True)

    for sub in sub_dirs:
        os.makedirs(os.path.join(main_dir, sub), exist_ok=True)


def fragment_genome(input_fasta, output_fasta, fragment_size=1020):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            length = len(seq)

            for i in range(0, length, fragment_size):
                fragment = seq[i:i+fragment_size]

                if len(fragment) == fragment_size:
                    header = f"{record.id}_frag_{i}"
                    out.write(f">{header}\n{fragment}\n")


def make_blast_db(fasta, db_path):
    cmd = [
        "makeblastdb",
        "-in", fasta,
        "-dbtype", "nucl",
        "-out", db_path
    ]
    subprocess.run(cmd, check=True)


def run_blast(query, db, output, threads):
    cmd = [
        "blastn",
        "-query", query,
        "-db", db,
        "-task", "blastn",
        "-dust", "no",
        "-xdrop_gap", "150",
        "-penalty", "-1",
        "-reward", "1",
        "-evalue", "1e-15",
        "-max_target_seqs", "1",
        "-outfmt", "6",
        "-num_threads", str(threads),
        "-out", output
    ]

    subprocess.run(cmd, check=True)


def calculate_directional_ani(blast_file, fragment_size=1020):

    headers = ["query", "subject", "pident", "alength", "mismatch",
               "gapopen", "qstart", "qend", "sstart", "send",
               "evalue", "bitscore"]

    df = pd.read_csv(blast_file, sep="\t", header=None)
    df.columns = headers

    min_len = int(fragment_size * 0.35)
    df = df[df["alength"] >= min_len]

    df = df.loc[df.groupby("query")["bitscore"].idxmax()]

    if len(df) == 0:
        return 0

    return df["pident"].mean()


def calculate_reciprocal_ani(fwd_file, rev_file, fragment_size=1020):

    headers = ["query", "subject", "pident", "alength", "mismatch",
               "gapopen", "qstart", "qend", "sstart", "send",
               "evalue", "bitscore"]

    fwd = pd.read_csv(fwd_file, sep="\t", header=None)
    rev = pd.read_csv(rev_file, sep="\t", header=None)

    fwd.columns = headers
    rev.columns = headers

    min_len = int(fragment_size * 0.35)

    fwd = fwd[fwd["alength"] >= min_len]
    rev = rev[rev["alength"] >= min_len]

    fwd = fwd.loc[fwd.groupby("query")["bitscore"].idxmax()]
    rev = rev.loc[rev.groupby("query")["bitscore"].idxmax()]

    rbh = pd.merge(
        fwd,
        rev,
        left_on="subject",
        right_on="query",
        suffixes=("_fwd", "_rev")
    )

    rbh = rbh[
        (rbh["query_fwd"] == rbh["subject_rev"]) &
        (rbh["query_rev"] == rbh["subject_fwd"])
    ]

    if len(rbh) == 0:
        return 0

    rbh["identity_mean"] = (
        rbh["pident_fwd"] + rbh["pident_rev"]
    ) / 2

    return rbh["identity_mean"].mean()


# --------------------------------------------------
# Main Pipeline
# --------------------------------------------------

def main():

    args = parse_args()

    genomes = get_genome_files(args.input)

    if len(genomes) < 2:
        raise ValueError("Need at least 2 genome files.")

    results = []

    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):

            g1 = genomes[i]
            g2 = genomes[j]

            name1 = os.path.splitext(os.path.basename(g1))[0]
            name2 = os.path.splitext(os.path.basename(g2))[0]

            pair_dir = os.path.join(args.output, f"{name1}_vs_{name2}")
            create_structure(pair_dir)

            fragA = os.path.join(pair_dir, "fragment/A_frag.fasta")
            fragB = os.path.join(pair_dir, "fragment/B_frag.fasta")

            fragment_genome(g1, fragA)
            fragment_genome(g2, fragB)

            dbA = os.path.join(pair_dir, "mkdb/A_db")
            dbB = os.path.join(pair_dir, "mkdb/B_db")

            make_blast_db(fragA, dbA)
            make_blast_db(fragB, dbB)

            fwd_out = os.path.join(pair_dir, "blast/A_vs_B.tsv")
            rev_out = os.path.join(pair_dir, "blast/B_vs_A.tsv")

            run_blast(fragA, dbB, fwd_out, args.threads)
            run_blast(fragB, dbA, rev_out, args.threads)

            ani_1_to_2 = calculate_directional_ani(fwd_out)
            ani_2_to_1 = calculate_directional_ani(rev_out)
            reciprocal_ani = calculate_reciprocal_ani(fwd_out, rev_out)

            results.append([
                name1,
                name2,
                ani_1_to_2,
                ani_2_to_1,
                reciprocal_ani
            ])

            print(f"[✔] Completed: {name1} vs {name2}")

    df = pd.DataFrame(
        results,
        columns=["Genome1", "Genome2", "ANI_1_to_2", "ANI_2_to_1", "Reciprocal_ANI"]
    )

    os.makedirs(args.output, exist_ok=True)
    output_file = os.path.join(args.output, "OrthoANI_results.tsv")

    df.to_csv(output_file, sep="\t", index=False)

    print("\n🎯 All comparisons completed.")
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
