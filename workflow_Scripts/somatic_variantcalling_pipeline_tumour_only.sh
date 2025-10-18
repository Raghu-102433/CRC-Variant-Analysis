#!/bin/bash

# Usage: ./tumor_only_vc_pipeline.sh <tumor_id>
# Example: ./tumor_only_vc_pipeline.sh SRR8993342

# ----------------------------
# INPUTS
# ----------------------------
TUMOR=$1

# ----------------------------
# PATHS & RESOURCES
# ----------------------------
WORKDIR="/media/raghu/new_volume/bioinfo/crc_variant_analysis/vc_work"
REF="$WORKDIR/ref_genome/Homo_sapiens_assembly38.fasta"
DBSNP="$WORKDIR/human_db/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INDELS="$WORKDIR/human_db/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS="$WORKDIR/human_db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
FUNCOTATOR="$WORKDIR/human_db/funcotator_dataSources.v1.8.hg38.20230908s"

cd "$WORKDIR"

# ----------------------------
# FASTQ FILES
# ----------------------------
FASTQ_TUMOR1="${TUMOR}_1.fastq.gz"
FASTQ_TUMOR2="${TUMOR}_2.fastq.gz"

# ----------------------------
# OUTPUT DIRECTORIES
# ----------------------------
OUTDIR="$WORKDIR/outputs/$TUMOR"
mkdir -p $OUTDIR/{qc_reports,trimmed_reads,bam_files,vcf_files,logs}

# ----------------------------
# STEP 1: QC
# ----------------------------
echo "Running FastQC..."
fastqc $FASTQ_TUMOR1 $FASTQ_TUMOR2 -o $OUTDIR/qc_reports

# ----------------------------
# STEP 2: Trimming
# ----------------------------
echo "Trimming..."
trimmomatic PE -threads 4 \
  $FASTQ_TUMOR1 $FASTQ_TUMOR2 \
  $OUTDIR/trimmed_reads/${TUMOR}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R1_unpaired.fq.gz \
  $OUTDIR/trimmed_reads/${TUMOR}_R2_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R2_unpaired.fq.gz \
  ILLUMINACLIP:$WORKDIR/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

fastqc -t 8 $OUTDIR/trimmed_reads/${TUMOR}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R2_paired.fq.gz \
  -o $OUTDIR/qc_reports

# ----------------------------
# STEP 3: Alignment
# ----------------------------
echo "Aligning..."
bwa mem -t 8 -R "@RG\tID:${TUMOR}\tSM:${TUMOR}\tPL:ILLUMINA" $REF \
  $OUTDIR/trimmed_reads/${TUMOR}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R2_paired.fq.gz \
  | samtools sort -o $OUTDIR/bam_files/${TUMOR}_sorted.bam

samtools index $OUTDIR/bam_files/${TUMOR}_sorted.bam

# ----------------------------
# STEP 4: Mark Duplicates
# ----------------------------
echo "Marking duplicates..."
gatk MarkDuplicatesSpark \
  -I $OUTDIR/bam_files/${TUMOR}_sorted.bam \
  -O $OUTDIR/bam_files/${TUMOR}_dedup.bam \
  -M $OUTDIR/logs/${TUMOR}_metrics.txt

# ----------------------------
# STEP 5: BQSR
# ----------------------------
echo "Base Quality Score Recalibration..."
gatk BaseRecalibrator \
  -I $OUTDIR/bam_files/${TUMOR}_dedup.bam -R $REF \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $MILLS \
  -O $OUTDIR/logs/${TUMOR}_recal.table

gatk ApplyBQSR \
  -I $OUTDIR/bam_files/${TUMOR}_dedup.bam -R $REF \
  --bqsr-recal-file $OUTDIR/logs/${TUMOR}_recal.table \
  -O $OUTDIR/bam_files/${TUMOR}_recal.bam

samtools index $OUTDIR/bam_files/${TUMOR}_recal.bam

# ----------------------------
# STEP 6: Somatic Variant Calling (Tumor-only)
# ----------------------------
RAW_VCF="$OUTDIR/vcf_files/${TUMOR}_tumor_only_raw.vcf.gz"
FILTERED_VCF="$OUTDIR/vcf_files/${TUMOR}_tumor_only_filtered.vcf.gz"
ANNOTATED_VCF="$OUTDIR/vcf_files/${TUMOR}_tumor_only_annotated.vcf.gz"

echo "Calling somatic variants with Mutect2 (tumor-only mode)..."
gatk --java-options "-Xmx16g" Mutect2 \
  -R $REF \
  -I $OUTDIR/bam_files/${TUMOR}_recal.bam -tumor $TUMOR \
  -O $RAW_VCF

# ----------------------------
# STEP 7: Filter Mutect2 Calls
# ----------------------------
echo "Filtering variants..."
gatk FilterMutectCalls \
  -R $REF \
  -V $RAW_VCF \
  -O $FILTERED_VCF

# Index filtered VCF
tabix -p vcf $FILTERED_VCF

# ----------------------------
# STEP 8: Annotation
# ----------------------------
echo "Annotating variants with Funcotator..."
gatk Funcotator \
  --variant $FILTERED_VCF \
  --reference $REF \
  --ref-version hg38 \
  --data-sources-path $FUNCOTATOR \
  --output $ANNOTATED_VCF \
  --output-file-format VCF

echo "Tumor-only Somatic Variant Calling Pipeline completed successfully!"
