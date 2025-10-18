#!/bin/bash

WORKDIR="/media/raghu/new_volume/bioinfo/crc_variant_analysis/vc_work/PON"
REF_GENOME_DIR="$WORKDIR/../ref_genome"
DB_DIR="$WORKDIR/../human_db"
REF="$REF_GENOME_DIR/Homo_sapiens_assembly38.fasta"
INTERVAL_FILE="$REF_GENOME_DIR/wgs_calling_regions.hg38.interval_list"
OUTDIR="$WORKDIR/PON_out"
mkdir -p "$OUTDIR"/{qc_reports,trimmed_reads,bam_files,vcf_files,logs}

NORMALS=("SRR8993331" "SRR8993344" "SRR8993385")

# =========================================================
# STEP 1: Prepare BAMs (trim, align, mark duplicates; no BQSR)
# =========================================================
echo "=== STEP 1: Preparing BAMs without BQSR ==="
for SAMPLE in "${NORMALS[@]}"; do
  echo "▶ Processing sample: $SAMPLE"

  FQ1="${SAMPLE}_1.fastq.gz"
  FQ2="${SAMPLE}_2.fastq.gz"
  BAM="$OUTDIR/bam_files/${SAMPLE}_dedup.bam"

  if [ -f "$BAM" ]; then
    echo "  → BAM exists, skipping."
    continue
  fi

  fastqc "$FQ1" "$FQ2" -o "$OUTDIR/qc_reports"

  trimmomatic PE -threads 4 \
    "$FQ1" "$FQ2" \
    "$OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz" "$OUTDIR/trimmed_reads/${SAMPLE}_R1_unpaired.fq.gz" \
    "$OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz" "$OUTDIR/trimmed_reads/${SAMPLE}_R2_unpaired.fq.gz" \
    ILLUMINACLIP:"$WORKDIR/../TruSeq3-PE.fa":2:30:10 LEADING:3 TRAILING:3 MINLEN:36

  bwa mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "$REF" \
    "$OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz" "$OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz" \
    | samtools sort -o "$OUTDIR/bam_files/${SAMPLE}_sorted.bam"

  samtools index "$OUTDIR/bam_files/${SAMPLE}_sorted.bam"

  gatk MarkDuplicatesSpark \
    -I "$OUTDIR/bam_files/${SAMPLE}_sorted.bam" \
    -O "$BAM" \
    -M "$OUTDIR/logs/${SAMPLE}_metrics.txt"

  samtools index "$BAM"
done

echo "✅ STEP 1 COMPLETE: BAMs ready."

# =========================================================
# STEP 2: Mutect2 Artifact Mode
# =========================================================
echo "=== STEP 2: Running Mutect2 artifact-only mode ==="
for SAMPLE in "${NORMALS[@]}"; do
  BAM="$OUTDIR/bam_files/${SAMPLE}_dedup.bam"
  ART_VCF="$OUTDIR/vcf_files/${SAMPLE}_artifacts.vcf.gz"

  if [ -f "$ART_VCF" ]; then
    echo "  → Artifact VCF exists, skipping."
    continue
  fi

  echo "  → Calling artifacts with Mutect2"
  gatk --java-options "-Xmx16g" Mutect2 \
    -R "$REF" \
    -I "$BAM" \
    -tumor "$SAMPLE" \
    --max-mnp-distance 0 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -O "$ART_VCF"
done

echo "✅ STEP 2 COMPLETE: Artifact VCFs ready."

# =========================================================
# STEP 3: Combine Artifact VCFs → Create PoN
# =========================================================
echo "=== STEP 3: Creating Panel of Normals (PoN) ==="

VCF_DIR="$OUTDIR/vcf_files"
DB_PATH="$OUTDIR/pon_genomicsdb"

VCF_ARGS=()
for SAMPLE in "${NORMALS[@]}"; do
  VCF_FILE="$VCF_DIR/${SAMPLE}_artifacts.vcf.gz"
  if [ ! -f "$VCF_FILE" ]; then
    echo "❌ Missing artifact VCF: $VCF_FILE"
    exit 1
  fi
  VCF_ARGS+=( -V "$VCF_FILE" )
done

# Import into GenomicsDB
gatk --java-options "-Xmx12g" GenomicsDBImport \
  -R "$REF" \
  -L "$INTERVAL_FILE" \
  --genomicsdb-workspace-path "$DB_PATH" \
  "${VCF_ARGS[@]}"

# Create PoN
gatk --java-options "-Xmx12g" CreateSomaticPanelOfNormals \
  -R "$REF" \
  -V "gendb://$DB_PATH" \
  -O "$VCF_DIR/PON.vcf.gz"

echo "✅ Panel of Normals successfully created at $VCF_DIR/PON.vcf.gz"

