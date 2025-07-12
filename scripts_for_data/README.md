# VEP Annotation Pipeline

This repository contains scripts to set up and run the Ensembl Variant Effect Predictor (VEP) with multiple plugins for annotating ClinVar variants using the GRCh38 reference genome.

## Repository Structure

| File                        | Description                                                             |
|-----------------------------|-------------------------------------------------------------------------|
| `vep_plugin_pipeline.sh`    | Downloads and configures VEP, cache, and required plugins               |
| `prepare_clinvar.sh`        | Prepares the ClinVar VCF file (e.g., filtering, compression)            |
| `vep_annotation_pipeline.sh`| Runs VEP annotation with the configured plugins                         |
| `filter_variants.sh`        | Filters variants based on clinical significance, review status, and more |
| `README.md`                 | Documentation and usage instructions                                    |

## Prerequisites

- macOS or Linux
- Docker (must be installed and running)
- `wget`, `bgzip`, `tabix` (can be installed via `htslib` or `samtools`)
- Internet access to download plugin and cache data
- Manual download access for some plugin datasets

## Setup Instructions

### Step 1: Install Docker

Please install Docker manually from the official website:
[Docker Desktop](https://www.docker.com/products/docker-desktop)

Once installed, start Docker Desktop before proceeding.

### Step 2: Prepare VEP Plugins and Cache

Run the following script to:
- Pull the Ensembl VEP Docker image
- Download the VEP cache (GRCh38)
- Download and process the following plugins:
  - REVEL
  - BayesDel (GRCh38 compatible)
  - VARITY
  - AlphaMissense
  - PrimateAI

```bash
./vep_plugin_pipeline.sh
```

**Note:** Some plugins require manual downloads. Ensure the required files are placed in the appropriate plugin directories as defined in the script.

### Step 3: Prepare ClinVar VCF

Run the following script to process and compress the ClinVar VCF file for use with VEP.

```bash
./prepare_clinvar.sh
```

Make sure to modify the script with the correct input path to your ClinVar file.

### Step 4: Annotate with VEP

Run the final annotation script that invokes the Docker container with all configured plugins.

```bash
./vep_annotation_pipeline.sh
```

## Plugin Notes

- **REVEL**: Converted CSV to TSV, sorted, bgzipped, and tabix-indexed
- **BayesDel**: GRCh38 sorted version used and renamed to .bed.gz with proper indexing
- **VARITY**: Extracted from tarball, sorted, bgzipped, and indexed
- **AlphaMissense**: Downloaded and indexed directly
- **PrimateAI**: Reformatted header, sorted numerically, bgzipped, and indexed

## Output

- The final annotated VCF will be output in compressed .vcf.gz format
- Intermediate files (cache, plugin data) are stored in
  ```bash
   ~/.vep/
  ```
   and plugin folders

## Additional Filtering

For additional filtering of variants based on clinical significance, review status, and other criteria, you can run the following script:
```bash
./filter_variants.sh
```
This script will filter the variants and generate a summary CSV file.

## Conclusion

This pipeline provides a comprehensive approach to annotating ClinVar variants using VEP and various plugins. Ensure that all prerequisites are met and follow the setup instructions carefully for successful execution. If you encounter any issues, please refer to the documentation for troubleshooting tips or reach out for assistance.
