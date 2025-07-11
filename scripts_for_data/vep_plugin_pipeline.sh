#!/bin/bash

# Exit immediately if a command fails
set -e

# ----------------------------
# 1. INSTALL DOCKER (Manual Step for Mac)
# ----------------------------
echo "Please install Docker manually from https://www.docker.com/products/docker-desktop"
echo "Then launch Docker Desktop before continuing."
sleep 5

# ----------------------------
# 2. PULL VEP DOCKER IMAGE
# ----------------------------
docker pull ensemblorg/ensembl-vep

# ----------------------------
# 3. CREATE VEP DIRECTORY STRUCTURE
# ----------------------------
mkdir -p ~/.vep/Plugins
mkdir -p ~/.vep/homo_sapiens/114_GRCh38
cd ~/.vep/homo_sapiens/114_GRCh38

# ----------------------------
# 3A. DOWNLOAD CACHE FILES USING DOCKER
# ----------------------------
docker run -t -i -v ~/.vep:/data ensemblorg/ensembl-vep \
    vep_install -a cf -s homo_sapiens -y GRCh38 -c /data --assembly GRCh38

# ----------------------------
# 4. DOWNLOAD REFERENCE FASTA
# ----------------------------
cd ~/.vep/homo_sapiens/114_GRCh38
wget ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai

# ----------------------------
# 5. CLONE OFFICIAL VEP PLUGINS
# ----------------------------
cd ~/.vep/Plugins
git clone https://github.com/Ensembl/VEP_plugins.git .

# ----------------------------
# 6A. REVEL PLUGIN
# ----------------------------
wget https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "."' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz

# ----------------------------
# 6B. BAYESDEL PLUGIN (Revamp for GRCh38)
# ----------------------------
echo "==> Preparing BayesDel plugin"
PLUGIN_DIR="$HOME/.vep/Plugins"
mkdir -p "$PLUGIN_DIR/BayesDel_170824/BayesDel_170824_addAF"
# Ensure file is downloaded manually from your Google Drive
mv BayesDel_GRCh38_sorted.txt.gz "$PLUGIN_DIR/BayesDel_170824/BayesDel_170824_addAF/bayesdel_38.sorted.bed.gz"
tabix -f -p bed "$PLUGIN_DIR/BayesDel_170824/BayesDel_170824_addAF/bayesdel_38.sorted.bed.gz"

# ----------------------------
# 6C. VARITY PLUGIN
# ----------------------------
mkdir -p VARITY
cd VARITY
wget http://varity.varianteffect.org/downloads/varity_all_predictions.tar.gz
tar -xzvf varity_all_predictions.tar.gz
cat varity_all_predictions.txt | (head -n 1 && tail -n +2 | sort -t$'\t' -k 1,1 -k 2,2n) > varity_all_predictions_sorted.tsv
sed '1s/.*/#&/' varity_all_predictions_sorted.tsv > varity_all_predictions.tsv
bgzip varity_all_predictions.tsv
tabix -f -s 1 -b 2 -e 2 varity_all_predictions.tsv.gz
cd ..

# ----------------------------
# 6D. ALPHAMISSENSE PLUGIN
# ----------------------------
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz

# ----------------------------
# 6E. PRIMATEAI PLUGIN
# ----------------------------
gunzip -cf PrimateAI_scores_v0.2_hg38.tsv.gz | \
  sed '12s/.*/#&/' | sed '/^$/d' | \
  awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k2,2n -V"}' | \
  bgzip > PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz
tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz

# ----------------------------
# 7. COMPRESS INPUT VCF FILE
# ----------------------------
bgzip /Users/sharmishtaganesh/Desktop/SEM2/ML_Sem2Proj/SnpEff_Output/clinvar.ann.vcf
tabix -p vcf /Users/sharmishtaganesh/Desktop/SEM2/ML_Sem2Proj/SnpEff_Output/clinvar.ann.vcf.gz

echo "All plugin preparation steps completed successfully."
