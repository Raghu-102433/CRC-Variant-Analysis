#!/bin/bash
# Usage: ./run_vcf_analysis.sh <input_vep_vcf> <output_dir>
# Example: ./run_vcf_analysis.sh input.vep.vcf.gz analysis/

set -euo pipefail

VCF=$1
OUTDIR=$2
TABLE="${OUTDIR}/variants.table"

# Ensure output dir exists
mkdir -p "$OUTDIR"

# Step 1: Convert VCF to table
echo "Converting VEP VCF to table..."
gatk VariantsToTable \
    -V "$VCF" \
    -F CHROM -F POS -F REF -F ALT -F CSQ -F TLOD -F NLOD -F AF \
    -O "$TABLE"

echo "VCF converted to table at $TABLE"

# Step 2: Run Python plotting script
echo "Running variant analysis and plotting..."
python3 vcf_analyser.py "$TABLE" "$OUTDIR"

echo "✅ Analysis completed. Results in $OUTDIR"

