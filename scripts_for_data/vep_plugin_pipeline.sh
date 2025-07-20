#!/bin/bash

set -e  # Exit if a command fails

# Set base directory on external drive
BASE="/Volumes/Sharmis_drive"
VEP_DIR="$BASE/.vep"
PLUGIN_DIR="$VEP_DIR/Plugins"

# 1. Docker Check
echo "Please install Docker manually from https://www.docker.com/products/docker-desktop"
echo "Then launch Docker Desktop before continuing."
sleep 5

# 2. Pull Ensembl VEP Docker Image
docker pull ensemblorg/ensembl-vep

# 3. Setup VEP Cache & FASTA
mkdir -p "$VEP_DIR/homo_sapiens/114_GRCh38"
mkdir -p "$PLUGIN_DIR"
cd "$VEP_DIR/homo_sapiens/114_GRCh38"

docker run -t -i -v "$VEP_DIR:/data" ensemblorg/ensembl-vep \
    vep_install -a cf -s homo_sapiens -y GRCh38 -c /data --assembly GRCh38

# Reference FASTA
wget -c ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai

# 4. VEP Official Plugins
cd "$PLUGIN_DIR"
git clone https://github.com/Ensembl/VEP_plugins.git .

# 5. REVEL Plugin
cd "$PLUGIN_DIR"
wget https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "."' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz

# 6. BayesDel Plugin
mkdir -p "$PLUGIN_DIR/BayesDel_170824/BayesDel_170824_addAF"
mv BayesDel_GRCh38_sorted.txt.gz "$PLUGIN_DIR/BayesDel_170824/BayesDel_170824_addAF/bayesdel_38.sorted.bed.gz"
tabix -f -p bed "$PLUGIN_DIR/BayesDel_170824/BayesDel_170824_addAF/bayesdel_38.sorted.bed.gz"

# 7. VARITY Plugin
mkdir -p "$PLUGIN_DIR/VARITY"
cd "$PLUGIN_DIR/VARITY"
wget http://varity.varianteffect.org/downloads/varity_all_predictions.tar.gz
tar -xzvf varity_all_predictions.tar.gz
cat varity_all_predictions.txt | (head -n 1 && tail -n +2 | sort -t$'\t' -k 1,1 -k 2,2n) > varity_all_predictions_sorted.tsv
sed '1s/.*/#&/' varity_all_predictions_sorted.tsv > varity_all_predictions.tsv
bgzip varity_all_predictions.tsv
tabix -f -s 1 -b 2 -e 2 varity_all_predictions.tsv.gz

# 8. AlphaMissense Plugin
cd "$PLUGIN_DIR"
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz

# 9. PrimateAI Plugin
gunzip -cf PrimateAI_scores_v0.2_hg38.tsv.gz | \
  sed '12s/.*/#&/' | sed '/^$/d' | \
  awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k2,2n -V"}' | \
  bgzip > PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz
tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz

# 10. dbscSNV Plugin
cd "$PLUGIN_DIR"
wget https://sites.google.com/site/jpopgen/dbNSFP/dbscSNV1.1.zip
unzip dbscSNV1.1.zip
bgzip dbscSNV1.1_GRCh38.txt
tabix -s 1 -b 2 -e 2 dbscSNV1.1_GRCh38.txt.gz

# 11. dbNSFP (GRCh38)
mkdir -p "$BASE/dbNSFP"
cd "$BASE/dbNSFP"
wget -c https://download.genos.us/dbnsfp/academic/ac0ddc/BGZF_format/dbNSFP5.1a_grch38.gz

# 12. KOVA Dataset
KOVA_DIR="$BASE/KOVA"
mkdir -p "$KOVA_DIR"
for chr in {1..22} X Y; do
    url="https://www.kobic.re.kr/kova/file_downloads?name=1_KOVA.v7.chr${chr}.tsv.gz"
    out="${KOVA_DIR}/KOVA.v7.chr${chr}.tsv.gz"
    echo "Downloading chr${chr}..."
    wget -c -O "$out" "$url"
done

# 13. dbSNP b151 GRCh38.p13
mkdir -p "$BASE/dbSNP"
cd "$BASE/dbSNP"
wget -c -O dbSNP_b151_GRCh38.vcf.gz \
  https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
wget -c -O dbSNP_b151_GRCh38.vcf.gz.tbi \
  https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi

# 14. gnomAD Exomes
GNOMAD_DIR="$BASE/gnomAD"
mkdir -p "$GNOMAD_DIR"
cd "$GNOMAD_DIR"
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi


# All Done

echo "All VEP plugins, data sources, and input VCF file setup completed under $BASE."
