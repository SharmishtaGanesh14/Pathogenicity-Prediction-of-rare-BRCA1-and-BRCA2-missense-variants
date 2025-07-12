#!/bin/bash


# Variant Annotation Pipeline using VEP (Docker)
# Author: Sharmishta
# Date: July 2025
# Purpose: Annotate ClinVar VCF using Ensembl VEP via Docker


# Step 1: Set Variables

# Path to your VCF file
INPUT_VCF="clinvar.ann.vcf"
CHUNK_PREFIX="chunk"
MAX_CHUNK_SIZE=1000000000 # 1GB
VEP_OUTPUT_DIR="/project/ML_Sem2Proj/SnpEff_Output"
DOCKER_VEP_IMG="ensemblorg/ensembl-vep"
CACHE_DIR="$HOME/.vep"

# Step 2: Prepare VEP Data (once) 
# These files are assumed to already be downloaded to ~/.vep
# - Cache: homo_sapiens/114_GRCh38
# - FASTA: Homo_sapiens.GRCh38.dna.toplevel.fa.gz
# - Plugins: AlphaMissense, BayesDel, PrimateAI, REVEL

# Step 3: Split VCF File into 1GB Chunks 

echo "Splitting large VCF into chunks of 1GB..."

# Extract header
grep '^#' "$INPUT_VCF" > header.vcf

# Split non-header part by byte size
grep -v '^#' "$INPUT_VCF" | split -b $MAX_CHUNK_SIZE -d -a 3 - "$CHUNK_PREFIX.data."

# Reattach header to each chunk and compress/index
for chunk in ${CHUNK_PREFIX}.data.*; do
    newfile="${CHUNK_PREFIX}_$(basename $chunk).vcf"
    cat header.vcf "$chunk" > "$newfile"

    # Compress and index
    bgzip -f "$newfile"
    tabix -p vcf "${newfile}.gz"
    rm "$chunk"
done

rm header.vcf
echo "VCF chunks created, compressed, and indexed."

# Step 4: Run VEP via Docker for Each Chunk

echo "Running VEP annotation on each VCF chunk..."

for gzvcf in ${CHUNK_PREFIX}_*.vcf.gz; do
    base=$(basename "$gzvcf" .vcf.gz)
    echo "▶ Annotating $gzvcf..."

    docker run -t -i \
      -v "$CACHE_DIR":/data \
      -v "$(pwd)":/project \
      $DOCKER_VEP_IMG \
      vep \
        -i /project/$gzvcf \
        -o /project/${base}.vep.vcf.gz \
        --cache --offline \
        --dir_cache /data \
        --dir_plugins /data/Plugins \
        --fasta /data/homo_sapiens/114_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --assembly GRCh38 \
        --format vcf --vcf \
        --compress_output bgzip \
        --force_overwrite \
        --symbol --canonical --protein --biotype --uniprot --hgvs \
        --no_stats \
        --plugin AlphaMissense,file=/data/Plugins/AlphaMissense_hg38.tsv.gz \
        --plugin BayesDel,file=/data/Plugins/BayesDel_170824/BayesDel_170824_addAF/bayesdel_38.sorted.bed.gz \
        --plugin PrimateAI,/data/Plugins/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz \
        --plugin REVEL,/data/Plugins/new_tabbed_revel.tsv.gz

    # Index output
    tabix -p vcf "/project/${base}.vep.vcf.gz"
done

echo "All chunks annotated and indexed."

# Optional: Merge annotated files later using bcftools 
# bcftools concat -a *.vep.vcf.gz -Oz -o merged_annotated.vcf.gz

# END OF PIPELINE

# Summary of Files and Their Purposes
# | File/Directory                                               | Purpose                                       |
# |--------------------------------------------------------------|-----------------------------------------------|
# | `clinvar.ann.vcf`                                            | Input VCF from SnpEff                         |
# | `~/.vep/homo_sapiens/114_GRCh38/`                            | VEP cache files                               |
# | `~/.vep/homo_sapiens/114_GRCh38/*.fa.gz`                     | Reference genome (FASTA) for VEP              |
# | `~/.vep/Plugins/`                                            | Directory containing plugin files             |
# | `AlphaMissense_hg38.tsv.gz`                                  | Plugin file for AlphaMissense                 |
# | `bayesdel_38.sorted.bed.gz`                                  | Plugin file for BayesDel (sorted BED)         |
# | `PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz`                | PrimateAI plugin file (bgzipped, tabix-indexed)|
# | `new_tabbed_revel.tsv.gz`                                    | REVEL plugin file (tabix-indexed TSV)         |

