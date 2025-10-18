#!/bin/bash

# Usage: ./somatic_vc_pipeline.sh <tumor_id> <normal_id>
# Example: ./somatic_vc_pipeline.sh SRR8993342 SRR8993344

set -euo pipefail

# Inputs
TUMOR=$1
NORMAL=$2

# Paths
WORKDIR="/media/raghu/new_volume/bioinfo/crc_variant_analysis/vc_work"
REF="$WORKDIR/ref_genome/Homo_sapiens_assembly38.fasta"
DBSNP="$WORKDIR/human_db/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INDELS="$WORKDIR/human_db/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS="$WORKDIR/human_db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
FUNCOTATOR="$WORKDIR/human_db/funcotator_dataSources.v1.8.hg38.20230908s

cd "$WORKDIR"

# FASTQ naming convention
FASTQ_TUMOR1="${TUMOR}_1.fastq.gz"
FASTQ_TUMOR2="${TUMOR}_2.fastq.gz"
FASTQ_NORMAL1="${NORMAL}_1.fastq.gz"
FASTQ_NORMAL2="${NORMAL}_2.fastq.gz"

# Output structure (everything in outputs/)
OUTDIR="$WORKDIR/outputs"
mkdir -p $OUTDIR/{qc_reports,trimmed_reads,bam_files,vcf_files,logs}

# ============== PIPELINE ==============

echo "Running FastQC..."
fastqc $FASTQ_TUMOR1 $FASTQ_TUMOR2 $FASTQ_NORMAL1 $FASTQ_NORMAL2 -o $OUTDIR/qc_reports

echo "Trimming..."
trimmomatic PE -threads 4 \
  $FASTQ_TUMOR1 $FASTQ_TUMOR2 \
  $OUTDIR/trimmed_reads/${TUMOR}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R1_unpaired.fq.gz \
  $OUTDIR/trimmed_reads/${TUMOR}_R2_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R2_unpaired.fq.gz \
  ILLUMINACLIP:$WORKDIR/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

trimmomatic PE -threads 4 \
  $FASTQ_NORMAL1 $FASTQ_NORMAL2 \
  $OUTDIR/trimmed_reads/${NORMAL}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${NORMAL}_R1_unpaired.fq.gz \
  $OUTDIR/trimmed_reads/${NORMAL}_R2_paired.fq.gz $OUTDIR/trimmed_reads/${NORMAL}_R2_unpaired.fq.gz \
  ILLUMINACLIP:$WORKDIR/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
 
echo "Running FastQC on trimmed reads..."
for fq in \
  $OUTDIR/trimmed_reads/${TUMOR}_R1_paired.fq.gz \
  $OUTDIR/trimmed_reads/${TUMOR}_R2_paired.fq.gz \
  $OUTDIR/trimmed_reads/${NORMAL}_R1_paired.fq.gz \
  $OUTDIR/trimmed_reads/${NORMAL}_R2_paired.fq.gz
do
  fastqc -t 8 "$fq" -o $OUTDIR/qc_reports
done

echo "Aligning..."
bwa mem -t 8 -R "@RG\tID:${TUMOR}\tSM:${TUMOR}\tPL:ILLUMINA" $REF \
  $OUTDIR/trimmed_reads/${TUMOR}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R2_paired.fq.gz \
  | samtools sort -o $OUTDIR/bam_files/${TUMOR}_sorted.bam

bwa mem -t 8 -R "@RG\tID:${NORMAL}\tSM:${NORMAL}\tPL:ILLUMINA" $REF \
  $OUTDIR/trimmed_reads/${NORMAL}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${NORMAL}_R2_paired.fq.gz \
  | samtools sort -o $OUTDIR/bam_files/${NORMAL}_sorted.bam

samtools index $OUTDIR/bam_files/${TUMOR}_sorted.bam
samtools index $OUTDIR/bam_files/${NORMAL}_sorted.bam

echo "Marking duplicates..."
gatk MarkDuplicatesSpark \
  -I $OUTDIR/bam_files/${TUMOR}_sorted.bam \
  -O $OUTDIR/bam_files/${TUMOR}_dedup.bam \
  -M $OUTDIR/logs/${TUMOR}_metrics.txt

gatk MarkDuplicatesSpark \
  -I $OUTDIR/bam_files/${NORMAL}_sorted.bam \
  -O $OUTDIR/bam_files/${NORMAL}_dedup.bam \
  -M $OUTDIR/logs/${NORMAL}_metrics.txt

echo "BQSR..."
gatk BaseRecalibrator \
  -I $OUTDIR/bam_files/${TUMOR}_dedup.bam -R $REF \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $MILLS \
  -O $OUTDIR/logs/${TUMOR}_recal.table
  
gatk ApplyBQSR \
  -I $OUTDIR/bam_files/${TUMOR}_dedup.bam -R $REF \
  --bqsr-recal-file $OUTDIR/logs/${TUMOR}_recal.table \
  -O $OUTDIR/bam_files/${TUMOR}_recal.bam
samtools index $OUTDIR/bam_files/${TUMOR}_recal.bam

gatk BaseRecalibrator \
  -I $OUTDIR/bam_files/${NORMAL}_dedup.bam -R $REF \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $MILLS \
  -O $OUTDIR/logs/${NORMAL}_recal.table
  
gatk ApplyBQSR \
  -I $OUTDIR/bam_files/${NORMAL}_dedup.bam -R $REF \
  --bqsr-recal-file $OUTDIR/logs/${NORMAL}_recal.table \
  -O $OUTDIR/bam_files/${NORMAL}_recal.bam
samtools index $OUTDIR/bam_files/${NORMAL}_recal.bam

RAW_VCF="$OUTDIR/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_raw.vcf.gz"
FILTERED_VCF="$OUTDIR/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_filtered.vcf.gz"
ANNOTATED_VCF="$OUTDIR/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_annotated.vcf.gz"

echo "Calling somatic variants..."
gatk --java-options "-Xmx12g" Mutect2 \
  -R $REF \
  -I $OUTDIR/bam_files/${TUMOR}_recal.bam -tumor $TUMOR \
  -I $OUTDIR/bam_files/${NORMAL}_recal.bam -normal $NORMAL \
  -O $OUTDIR/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_raw.vcf.gz

echo "Filtering Mutect2 calls..."
gatk FilterMutectCalls \
  -R $REF \
  -V $RAW_VCF \
  -O $FILTERED_VCF

# Index the filtered VCF
if [ ! -f "${FILTERED_VCF}.tbi" ]; then
    tabix -p vcf $FILTERED_VCF
fi

echo "Annotating with Funcotator..."
gatk Funcotator \
  --variant $FILTERED_VCF \
  --reference $REF \
  --ref-version hg38 \
  --data-sources-path $FUNCOTATOR \
  --output $ANNOTATED_VCF \
  --output-file-format VCF
  
echo "Somatic Variant Calling Pipeline completed successfully!" 
