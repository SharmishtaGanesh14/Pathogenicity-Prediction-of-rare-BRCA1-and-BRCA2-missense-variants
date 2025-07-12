#!/bin/bash
set -euo pipefail

# Configuration
INPUT_VCF="merged_annotated.vep.vcf.gz"
FINAL_VCF="filtered_final.vcf.gz"
CLIN_SIG_VALUES="Benign,Benign/Likely_benign,Likely_benign,Pathogenic,Pathogenic/Likely_pathogenic,Likely_pathogenic"
REVIEW_STATUS="practice_guideline,reviewed_by_expert_panel,criteria_provided_multiple_submitters_no_conflicts"
CANCER_GENES="APC,ATM,BAP1,BARD1,BRCA1,BRCA2,BRIP1,CDH1,CDK4,CDKN2A,CHEK2,EPCAM,MITF,MLH1,MSH2,MSH6,MUTYH,NBN,PALB2,PMS2,POLD1,POLE,PTEN,RAD51C,RAD51D,SMAD4,STK11,TP53"

# Step 1: Filter by clinical significance
echo "STEP 1/6: Filtering by clinical significance..."
bcftools view -i "CLNSIG ~ '${CLIN_SIG_VALUES//,/\|}'" ${INPUT_VCF} -Oz -o step1_clinsig.vcf.gz
bcftools index step1_clinsig.vcf.gz

# Step 2: Filter by ClinVar review status  
echo "STEP 2/6: Filtering by review status..."
bcftools view -i "CLNREVSTAT ~ '${REVIEW_STATUS//,/\|}'" step1_clinsig.vcf.gz -Oz -o step2_reviewed.vcf.gz
bcftools index step2_reviewed.vcf.gz

# Step 3: Filter variant type (SNVs only)
echo "STEP 3/6: Filtering variant type (SNVs only)..."
bcftools view -i 'TYPE="snp"' step2_reviewed.vcf.gz -Oz -o step3_snvs.vcf.gz
bcftools index step3_snvs.vcf.gz

# Step 4: Filter by MAF (<0.005)
echo "STEP 4/6: Filtering by MAF (<0.005)..."
bcftools view -e 'INFO/AF >= 0.005' step3_snvs.vcf.gz -Oz -o step4_maf.vcf.gz
bcftools index step4_maf.vcf.gz

# Step 5: Filter by consequence (missense/splice)
echo "STEP 5/6: Filtering by consequence..."
bcftools view -i 'CSQ ~ "missense_variant" || CSQ ~ "missense_variant&splice_region_variant"' step4_maf.vcf.gz -Oz -o step5_missense.vcf.gz
bcftools index step5_missense.vcf.gz

# Step 6: Filter hereditary cancer genes
echo "STEP 6/6: Filtering cancer genes..."
bcftools view -i "GENEINFO ~ '${CANCER_GENES//,/\|}'" step5_missense.vcf.gz -Oz -o ${FINAL_VCF}
bcftools index ${FINAL_VCF}

# Generate summary statistics
echo "Generating summary..."
bcftools query -f '%INFO/GENEINFO\t%INFO/CLNSIG\n' ${FINAL_VCF} | \
awk -F '\t' '
BEGIN {
    OFS=","
    print "Gene,P/LP,B/LB,All"
} 
{
    split($1,gene,":")
    clnsig=$2
    
    # Count classifications
    p_lp = (clnsig~/Pathogenic/||clnsig~/Likely_pathogenic/) ? 1 : 0
    b_lb = (clnsig~/Benign/||clnsig~/Likely_benign/) ? 1 : 0

    counts[gene[1]]++
    if(p_lp) p_lp_counts[gene[1]]++
    if(b_lb) b_lb_counts[gene[1]]++
}
END {
    for(gene in counts) {
        print gene, p_lp_counts[gene]+0, b_lb_counts[gene]+0, counts[gene]
    }
}' | sort > "supplementary_materials/gene_clnsig_summary.csv"

echo "Pipeline completed successfully!"
echo "Final filtered VCF: ${FINAL_VCF}"
echo "Summary CSV created: supplementary_materials/gene_clnsig_summary.csv"

# Clean up intermediate files
rm step*.vcf.gz*
