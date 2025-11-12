
#!/bin/bash

# SparCC Network Analysis Pipeline
# This script runs the complete SparCC analysis including correlation calculation and significance testing

echo "=== Starting SparCC Network Analysis Pipeline ==="

# Set up directories
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"

# Create necessary directories
mkdir -p bootstraps
mkdir -p pvals
mkdir -p logs
mkdir -p results

echo "1. Setting up environment..."
export PATH="$HOME/qyy/BPD1/Sparcc/SparCC3:$PATH"

# Check if SparCC is available
if ! command -v SparCC.py &> /dev/null; then
    echo "Error: SparCC.py not found in PATH"
    echo "Please ensure SparCC is installed and added to PATH"
    exit 1
fi

echo "2. Running SparCC correlation calculation..."

# Define groups to analyze
groups=("none" "mild" "moderate" "severe")

for group in "${groups[@]}"; do
    echo "Processing group: $group"
    
    # Create group-specific directories
    mkdir -p "$group"
    mkdir -p "bootstraps/$group"
    mkdir -p "pvals/$group"
    
    input_file="$group/norm_rarefy_${group}_relabundance.txt"
    
    if [[ ! -f "$input_file" ]]; then
        echo "Warning: Input file not found: $input_file"
        continue
    fi
    
    echo "  Input file: $input_file"
    
    # Step 1: Calculate correlation matrix
    echo "  Calculating correlation matrix..."
    SparCC.py "$input_file" \
        -i 20 \
        --cor_file="$group/cor_sparcc.norm_rarefy_${group}.txt" \
        >> "logs/sparcc_${group}.log" 2>&1
    
    # Step 2: Generate bootstrap datasets
    echo "  Generating bootstrap datasets..."
    MakeBootstraps.py "$input_file" \
        -n 100 \
        -t "bootstraps/${group}/${group}_bootstrap_#.txt" \
        -p "pvals/${group}/" \
        >> "logs/bootstrap_${group}.log" 2>&1
    
    # Step 3: Calculate correlations for bootstrap datasets
    echo "  Calculating bootstrap correlations..."
    for bootstrap_file in bootstraps/${group}/${group}_bootstrap_*.txt; do
        if [[ -f "$bootstrap_file" ]]; then
            base_name=$(basename "$bootstrap_file" .txt)
            SparCC.py "$bootstrap_file" \
                --cor_file="bootstraps/${group}/cor_${base_name}.txt" \
                -i 10 \
                >> "logs/bootstrap_cor_${group}.log" 2>&1
        fi
    done
    
    # Step 4: Calculate p-values
    echo "  Calculating p-values..."
    PseudoPvals.py "$group/cor_sparcc.norm_rarefy_${group}.txt" \
        "bootstraps/${group}/cor_${group}_bootstrap_*.txt" \
        100 \
        -o "pvals/${group}/pvals.${group}.txt" \
        -t two_sided \
        >> "logs/pvals_${group}.log" 2>&1
    
    echo "  Completed: $group"
done

echo "3. Generating summary report..."

# Create analysis summary
summary_file="results/sparcc_analysis_summary.txt"
echo "SparCC Analysis Summary" > "$summary_file"
echo "Generated: $(date)" >> "$summary_file"
echo "=========================" >> "$summary_file"

for group in "${groups[@]}"; do
    echo "" >> "$summary_file"
    echo "Group: $group" >> "$summary_file"
    
    cor_file="$group/cor_sparcc.norm_rarefy_${group}.txt"
    pval_file="pvals/${group}/pvals.${group}.txt"
    
    if [[ -f "$cor_file" ]]; then
        echo "  Correlation matrix: ✓ ($(wc -l < "$cor_file") ASVs)" >> "$summary_file"
    else
        echo "  Correlation matrix: ✗" >> "$summary_file"
    fi
    
    if [[ -f "$pval_file" ]]; then
        echo "  P-values: ✓" >> "$summary_file"
    else
        echo "  P-values: ✗" >> "$summary_file"
    fi
done

echo "" >> "$summary_file"
echo "Output directories:" >> "$summary_file"
echo "  - bootstraps/    : Bootstrap datasets and correlations" >> "$summary_file"
echo "  - pvals/         : P-value results" >> "$summary_file"
echo "  - logs/          : Analysis logs" >> "$summary_file"
echo "  - results/       : Summary reports" >> "$summary_file"

echo "=== SparCC Analysis Pipeline Completed ==="
echo "Summary saved to: $summary_file"
