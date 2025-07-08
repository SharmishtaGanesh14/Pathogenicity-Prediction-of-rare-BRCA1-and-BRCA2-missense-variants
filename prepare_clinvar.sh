#!/bin/bash

# Author: Sharmishta Ganesh
# Date: 2025-07-08
# Description: This script downloads the latest ClinVar VCF for GRCh38 (RefSeq),
#              processes the GRCh38 reference genome to remove unwanted 'chr' prefixes
#              and ensure contig names match those used in ClinVar VCF, and finally
#              normalizes the variants using bcftools.

set -e  # Exit immediately if a command exits with a non-zero status.

# ----------------------------------------
# Step 1: Setup Working Directory
# ----------------------------------------
mkdir -p ClinVar_GRCh38 && cd ClinVar_GRCh38

# ----------------------------------------
# Step 2: Download ClinVar VCF (GRCh38)
# Source: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# ----------------------------------------
echo "Downloading ClinVar VCF..."
wget -N https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -N https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# ----------------------------------------
# Step 3: Download GRCh38.p14 Reference Genome (RefSeq build)
# Source: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
# ----------------------------------------
echo "Downloading GRCh38.p14 reference genome..."
wget -N https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# ----------------------------------------
# Step 4: Preprocess and Index the Reference Genome
# ----------------------------------------
echo "Unzipping and processing reference..."
gunzip -c GCF_000001405.40_GRCh38.p14_genomic.fna.gz > GRCh38_refseq_raw.fna

echo "Cleaning 'chr' and fixing MT naming..."
sed -E 's/^>chr([0-9XY]+)\s.*/>\1/; s/^>chrM\s.*/>MT/' GRCh38_refseq_raw.fna > GRCh38_no_chr.fna

echo "Indexing the cleaned reference..."
samtools faidx GRCh38_no_chr.fna

# ----------------------------------------
# Step 5: Normalize ClinVar VCF using bcftools
# ----------------------------------------
echo "Normalizing ClinVar VCF using bcftools..."
bcftools norm \
  -f GRCh38_no_chr.fna \
  clinvar.vcf.gz \
  -Oz -o clinvar.norm.vcf.gz

# ----------------------------------------
# Step 6: Index the Normalized VCF
# ----------------------------------------
echo "Indexing normalized VCF..."
tabix -p vcf clinvar.norm.vcf.gz

# ----------------------------------------
# Completion Message
# ----------------------------------------
echo "ClinVar VCF normalization complete!"
echo "Output VCF: clinvar.norm.vcf.gz"
echo "Reference used: GRCh38_no_chr.fna"
