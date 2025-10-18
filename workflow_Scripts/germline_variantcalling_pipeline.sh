#!/bin/bash

# Usage: ./germline_vc_pipeline.sh <sample_id>
# Example: ./germline_vc_pipeline.sh SRR8993342

set -euo pipefail

# Inputs
SAMPLE=$1

# Paths
WORKDIR="/media/raghu/new_volume/bioinfo/crc_variant_analysis/vc_work"
REF="$WORKDIR/ref_genome/Homo_sapiens_assembly38.fasta"
DBSNP="$WORKDIR/human_db/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INDELS="$WORKDIR/human_db/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS="$WORKDIR/human_db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
FUNCOTATOR="$WORKDIR/human_db/funcotator_dataSources.v1.8.hg38.20230908g"

cd "$WORKDIR"

# FASTQ naming convention
FASTQ1="${SAMPLE}_1.fastq.gz"
FASTQ2="${SAMPLE}_2.fastq.gz"

# Output structure
OUTDIR="$WORKDIR/outputs"
mkdir -p $OUTDIR/{qc_reports,trimmed_reads,bam_files,vcf_files,logs}

# ============== PIPELINE ==============

echo "Running FastQC..."
fastqc $FASTQ1 $FASTQ2 -o $OUTDIR/qc_reports

echo "Trimming..."
trimmomatic PE -threads 4 \
  $FASTQ1 $FASTQ2 \
  $OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${SAMPLE}_R1_unpaired.fq.gz \
  $OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz $OUTDIR/trimmed_reads/${SAMPLE}_R2_unpaired.fq.gz \
  ILLUMINACLIP:$WORKDIR/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

echo "Running FastQC on trimmed reads..."
for fq in \
  $OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz \
  $OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz
do
  fastqc -t 8 "$fq" -o $OUTDIR/qc_reports
done

echo "Aligning..."
bwa mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" $REF \
  $OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz \
  | samtools sort -o $OUTDIR/bam_files/${SAMPLE}_sorted.bam

samtools index $OUTDIR/bam_files/${SAMPLE}_sorted.bam

echo "Marking duplicates..."
gatk MarkDuplicatesSpark \
  -I $OUTDIR/bam_files/${SAMPLE}_sorted.bam \
  -O $OUTDIR/bam_files/${SAMPLE}_dedup.bam \
  -M $OUTDIR/logs/${SAMPLE}_metrics.txt

echo "BQSR..."
gatk BaseRecalibrator \
  -I $OUTDIR/bam_files/${SAMPLE}_dedup.bam -R $REF \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $MILLS \
  -O $OUTDIR/logs/${SAMPLE}_recal.table
  
gatk ApplyBQSR \
  -I $OUTDIR/bam_files/${SAMPLE}_dedup.bam -R $REF \
  --bqsr-recal-file $OUTDIR/logs/${SAMPLE}_recal.table \
  -O $OUTDIR/bam_files/${SAMPLE}_recal.bam

samtools index $OUTDIR/bam_files/${SAMPLE}_recal.bam

echo "Calling germline variants with HaplotypeCaller..."
RAW_VCF="$OUTDIR/vcf_files/${SAMPLE}_raw.vcf.gz"
gatk --java-options "-Xmx14g" HaplotypeCaller \
  -R $REF \
  -I $OUTDIR/bam_files/${SAMPLE}_recal.bam \
  -O $RAW_VCF \
  -ERC GVCF

# Paths for VCFs
RAW_VCF="$OUTDIR/vcf_files/${SAMPLE}_raw.vcf.gz"
GENOTYPED_VCF="$OUTDIR/vcf_files/${SAMPLE}_genotyped.vcf"
FILTERED_SNP_VCF="$OUTDIR/vcf_files/${SAMPLE}_filtered_snps.vcf"
FILTERED_INDEL_VCF="$OUTDIR/vcf_files/${SAMPLE}_filtered_indels.vcf"
MERGED_FILTERED_VCF="$OUTDIR/vcf_files/${SAMPLE}_filtered_merged.vcf.gz"
ANNOTATED_VCF="$OUTDIR/vcf_files/${SAMPLE}_annotated.vcf"

echo "Genotyping raw VCF..."
gatk GenotypeGVCFs \
  -R $REF \
  -V $RAW_VCF \
  -O $GENOTYPED_VCF

echo "Filtering SNPs..."
gatk SelectVariants \
  -R $REF \
  -V $GENOTYPED_VCF \
  --select-type-to-include SNP \
  -O "$OUTDIR/vcf_files/${SAMPLE}_raw_snps.vcf"

gatk VariantFiltration \
  -R $REF \
  -V "$OUTDIR/vcf_files/${SAMPLE}_raw_snps.vcf" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "SNP_Filter" \
  -O $FILTERED_SNP_VCF

echo "Filtering INDELs..."
gatk SelectVariants \
  -R $REF \
  -V $GENOTYPED_VCF \
  --select-type-to-include INDEL \
  -O "$OUTDIR/vcf_files/${SAMPLE}_raw_indels.vcf"

gatk VariantFiltration \
  -R $REF \
  -V "$OUTDIR/vcf_files/${SAMPLE}_raw_indels.vcf" \
  --filter-expression "QD < 2.0 || FS > 200.0" \
  --filter-name "INDEL_Filter" \
  -O $FILTERED_INDEL_VCF

echo "Merging filtered SNPs and INDELs..."
gatk MergeVcfs \
  -I $FILTERED_SNP_VCF \
  -I $FILTERED_INDEL_VCF \
  -O $MERGED_FILTERED_VCF

# Index merged VCF
if [ ! -f "${MERGED_FILTERED_VCF}.tbi" ]; then
    tabix -p vcf $MERGED_FILTERED_VCF
fi

echo "Annotating filtered variants with Funcotator..."
gatk Funcotator \
  --variant $MERGED_FILTERED_VCF \
  --reference $REF \
  --ref-version hg38 \
  --data-sources-path $FUNCOTATOR \
  --output $ANNOTATED_VCF \
  --output-file-format VCF
