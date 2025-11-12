#!/bin/bash
# SparCC correlation analysis workflow for BPD samples
# Date: 2025-11
# Requires: Python + SparCC + numpy

set -e
workdir=$(pwd)
input="otu_table.txt"

echo ">>> Running SparCC on ${input}"

# Step 1: Run SparCC
python SparCC.py $input --cor_file=sparcc.out

# Step 2: Generate bootstraps
mkdir -p bootstraps pvals
for i in {1..100}
do
    python MakeBootstraps.py $input -n 1 -t bootstraps/boot_$i.txt
done

# Step 3: Run SparCC on each bootstrap
for f in bootstraps/boot_*.txt
do
    python SparCC.py $f --cor_file=${f/.txt/.out}
done

# Step 4: Compute pseudo p-values
python PseudoPvals.py sparcc.out bootstraps/boot_*.out 100 -o pvals/pvals_two_sided.txt

# Step 5: Save log
date > sparcc.log
echo "SparCC pipeline finished successfully" >> sparcc.log
