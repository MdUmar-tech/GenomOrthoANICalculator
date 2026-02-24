# 🌐 GenomOrthoANICalculator – E2E Cloud Setup Guide

This guide explains how to:

- Connect to E2E cloud server
- Install dependencies
- Setup Python / Conda environment
- Run ANI pipeline safely
- Download results
- Clean server

---

# 🔹 1️⃣ Connect to Server

## Basic SSH Login

```bash
ssh -i id_rsa root@YOUR_SERVER_IP
```

Example:

```bash
ssh -i id_rsa root@164.52.219.234
```

## Debug SSH (if connection fails)

```bash
ssh -vvv root@164.52.209.136
```

---

# 🔹 2️⃣ Check Server Resources

Always check available resources before running heavy jobs.

```bash
nproc        # Number of CPU cores
lscpu        # Detailed CPU info
top          # Monitor CPU & RAM
df -h        # Disk usage
```

---

# 🐍 Option A: Setup Using Python Virtual Environment (venv)

## Update System

```bash
apt update
apt install python3 python3-pip -y
apt install ncbi-blast+ -y
```

## Create Virtual Environment

```bash
python3 -m venv bioenv
source bioenv/bin/activate
```

If you need to remove it:

```bash
rm -rf bioenv
```

## Install Required Packages

```bash
pip install pandas biopython
```

---

# 🧬 Option B: Setup Using Conda (Recommended)

## Install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

Check installation:

```bash
conda --version
```

---

## Create Environment

```bash
conda create -n ani_env -y
conda activate ani_env
```

Install required tools:

```bash
conda install -c bioconda ncbi-blast -y
conda install -c conda-forge pandas biopython -y
```

Verify BLAST:

```bash
blastn -version
```

---

# 📂 3️⃣ Upload Files to Server

From your local machine:

```bash
scp -i mykey.pem -r genomes GenomOrthoANICalculator.py root@164.52.219.234:/root/
```

---

# 🖥️ 4️⃣ Run Inside Screen (VERY IMPORTANT)

Start a screen session:

```bash
screen -S ani
```

Run the pipeline:

```bash
python GenomOrthoANICalculator.py \
--input genomes \
--output ANI_output \
--threads 20
```

---

## 🔹 Detach Screen (Keep Running)

Press:

```
Ctrl + A
then D
```

Reattach later:

```bash
screen -r ani
```

---

# 📥 5️⃣ Download Results

From your local machine:

```bash
scp -i mykey.pem -r root@164.52.219.234:/root/ANI_output/ .
```

---

# 🧹 6️⃣ Clean Server (Optional)

⚠ Only do this if fully finished.

```bash
rm -rf /root/*
```

Exit server:

```bash
exit
```

---

# ⚡ Recommended Thread Usage

Check CPU cores:

```bash
nproc
```

If server has 24 cores:

Use:

```bash
--threads 20
```

Always leave 2–4 cores free for stability.

---

# 🏆 Example Full Workflow

```bash
ssh -i id_rsa root@164.52.219.234
conda activate ani_env
screen -S ani
python GenomOrthoANICalculator.py --input genomes --output ANI_output --threads 20
```

---

# ✅ Best Practices

- ✔ Always use screen or tmux
- ✔ Check disk space before running
- ✔ Avoid spaces in genome filenames
- ✔ Download results before deleting data
- ✔ Monitor resources using top
- ✔ Use appropriate thread count

---

# 🎯 Summary

This guide ensures:

- Reproducible ANI analysis
- Proper cloud usage
- Efficient CPU utilization
- Safe long-running job execution
- Clean result retrieval

---

© Md Umar – GenomOrthoANICalculator