#!/bin/bash

# Author: Sharmishta Ganesh
# Date: 2025-07-08
# Description: This script downloads the latest ClinVar VCF for GRCh38 (RefSeq),
#              processes the GRCh38 reference genome (hs38d1 full + patch),
#              removes unwanted 'chr' prefixes, ensures contig name compatibility,
#              and finally normalizes the variants using bcftools.

set -e  # Exit on any command failure

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
# Step 3: Download GRCh38 Reference Genome (RefSeq + hs38d1 patches)
# Source: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/
# ----------------------------------------
echo "Downloading GRCh38 full + hs38d1 reference..."
wget -N https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz

# ----------------------------------------
# Step 4: Preprocess and Index the Reference Genome
# ----------------------------------------
echo "Processing reference genome..."
gunzip -c GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz > GRCh38_refseq_raw.fna

echo "Removing 'chr' prefix and fixing MT naming..."
sed -E 's/^>chr([0-9XY]+)\s.*/>\1/; s/^>chrM\s.*/>MT/' GRCh38_refseq_raw.fna > GRCh38_no_chr.fna

echo "[MANUAL STEP] Adding missing contigs (e.g., NT_113889.1, NW_009646201.1) to reference..."
echo "Ensure you've manually added missing contigs to GRCh38_no_chr.fna before proceeding."

# ----------------------------------------
# Step 5: Index the Cleaned Reference
# ----------------------------------------
echo "Indexing reference with samtools..."
samtools faidx GRCh38_no_chr.fna

# ----------------------------------------
# Step 6: Normalize ClinVar VCF using bcftools
# ----------------------------------------
echo "Normalizing ClinVar VCF using bcftools..."
bcftools norm \
  -f GRCh38_no_chr.fna \
  clinvar.vcf.gz \
  -Oz -o clinvar.norm.vcf.gz

# ----------------------------------------
# Step 7: Index the Normalized VCF
# ----------------------------------------
echo "Indexing normalized VCF..."
tabix -p vcf clinvar.norm.vcf.gz

# ----------------------------------------
# Completion Message
# ----------------------------------------
echo "ClinVar VCF normalization complete!"
echo "Output VCF: clinvar.norm.vcf.gz"
echo "Reference used: GRCh38_no_chr.fna (with manually added contigs)"
