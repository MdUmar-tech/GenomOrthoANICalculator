#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO

__version__ = "1.0.0"

# --------------------------------------------------
# Argument Parser
# --------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Optimized Pairwise OrthoANI calculation"
    )

    parser.add_argument("--input", required=True,
                        help="Folder containing genome fasta files")

    parser.add_argument("--output", required=True,
                        help="Output folder")

    parser.add_argument("--threads", type=int, default=2,
                        help="Threads for BLAST (default=2)")
   parser.add_argument( "--version", action="version", version="GenomOrthoANICalculator v1.0.0"
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


def create_main_structure(outdir):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "fragments"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "mkdb"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "blast"), exist_ok=True)


def fragment_genome(input_fasta, output_fasta, fragment_size=1020):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            length = len(seq)

            for i in range(0, length, fragment_size):
                fragment = seq[i:i + fragment_size]
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

    if os.path.getsize(blast_file) == 0:
        return 0

    headers = ["query", "subject", "pident", "alength", "mismatch",
               "gapopen", "qstart", "qend", "sstart", "send",
               "evalue", "bitscore"]

    df = pd.read_csv(blast_file, sep="\t", header=None)
    df.columns = headers

    min_len = int(fragment_size * 0.35)
    df = df[df["alength"] >= min_len]

    if len(df) == 0:
        return 0

    df = df.loc[df.groupby("query")["bitscore"].idxmax()]
    return df["pident"].mean()


def calculate_reciprocal_ani(fwd_file, rev_file, fragment_size=1020):

    if os.path.getsize(fwd_file) == 0 or os.path.getsize(rev_file) == 0:
        return 0

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

    if len(fwd) == 0 or len(rev) == 0:
        return 0

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

    rbh["identity_mean"] = (rbh["pident_fwd"] + rbh["pident_rev"]) / 2
    return rbh["identity_mean"].mean()


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():

    args = parse_args()
    genomes = get_genome_files(args.input)

    if len(genomes) < 2:
        raise ValueError("Need at least 2 genomes")

    create_main_structure(args.output)

    fragment_dir = os.path.join(args.output, "fragments")
    db_dir = os.path.join(args.output, "mkdb")
    blast_dir = os.path.join(args.output, "blast")

    print("🔹 Step 1: Fragmenting genomes & building BLAST DB (only once)")

    genome_map = {}

    for genome in genomes:
        name = os.path.splitext(os.path.basename(genome))[0]

        frag_file = os.path.join(fragment_dir, f"{name}_frag.fasta")
        db_file = os.path.join(db_dir, f"{name}_db")

        fragment_genome(genome, frag_file)
        make_blast_db(frag_file, db_file)

        genome_map[name] = {
            "fragment": frag_file,
            "db": db_file
        }

        print(f"[✔] Prepared: {name}")

    print("\n🔹 Step 2: Running pairwise BLAST")

    results = []
    names = list(genome_map.keys())

    for i in range(len(names)):
        for j in range(i + 1, len(names)):

            g1 = names[i]
            g2 = names[j]

            print(f"Comparing: {g1} vs {g2}")

            fwd_out = os.path.join(blast_dir, f"{g1}_vs_{g2}.tsv")
            rev_out = os.path.join(blast_dir, f"{g2}_vs_{g1}.tsv")

            run_blast(genome_map[g1]["fragment"],
                      genome_map[g2]["db"],
                      fwd_out,
                      args.threads)

            run_blast(genome_map[g2]["fragment"],
                      genome_map[g1]["db"],
                      rev_out,
                      args.threads)

            ani_1_to_2 = calculate_directional_ani(fwd_out)
            ani_2_to_1 = calculate_directional_ani(rev_out)
            reciprocal_ani = calculate_reciprocal_ani(fwd_out, rev_out)

            results.append([
                g1, g2,
                ani_1_to_2,
                ani_2_to_1,
                reciprocal_ani
            ])

    df = pd.DataFrame(results,
                      columns=["Genome1", "Genome2",
                               "ANI_1_to_2", "ANI_2_to_1",
                               "Reciprocal_ANI"])

    output_file = os.path.join(args.output, "OrthoANI_results.tsv")
    df.to_csv(output_file, sep="\t", index=False)

    print("\n🎯 All comparisons completed.")
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
