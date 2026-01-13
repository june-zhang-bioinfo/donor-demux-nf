# Donor Demultiplexing Nextflow Pipeline

A Nextflow pipeline for demultiplexing pooled single-cell RNA-seq samples using genetic variants (SNPs). The pipeline identifies genetically distinct donors in multiplexed 10x Genomics experiments through variant calling and probabilistic assignment, when no cell hashing (HTO) data is available. 

Genetically unrelated (or weakly related) donors and adequate SNP coverage present in both donor and pooled BAMs are required for the best performance.

## Overview

This pipeline performs:
1. **VCF Preparation** (optional): Processes donor BAM files to generate donor-specific variants
2. **Demultiplexing**: Uses Demuxalot to assign cells to donors based on genetic variants
3. **Results Visualization**: Generates PCA plots and sex-specific gene expression visualizations

## Pipeline Architecture

### Main Workflow (`main.nf`)
- Orchestrates the entire pipeline
- Supports two modes: full pipeline or work with pre-existing VCF

***Note on demo mode***
Donor genotyping BAMs are typically not publicly available due to privacy restrictions.
The demo mode allows users to run and validate the demultiplexing and QC components of the pipeline using a precomputed VCF.
The example VCF used for demo mode is derived from the Demuxafy documentation:
https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DataPrep.html

### Modules
- **`vcf_preparation.nf`**: Prepares donor variants from BAM files
  - `CHANGE_READ_GROUP`: Updates read group tags in BAM files
  - `FREEBAYES_CALL`: Performs variant calling using FreeBayes
- **`demultiplexing_analysis.nf`**: Performs demultiplexing and visualization
  - `DEMUXALOT`: Assigns cells to donors
  - `R_SEURAT_ANALYSIS`: Generates QC plots

## Requirements

### Software Dependencies
- **Singularity**
- **Nextflow** (24.04.2)

### Container Images
- `devjune26/donor_demux:1.0` - Contains samtools, bcftools, FreeBayes, and demuxalot
- `devjune26/r_eval:1.0` - Contains R with Seurat and visualization packages

## Pull Docker Images
```bash
mkdir containers
cd containers
module load singularity

# donor_demux container
singularity pull donor_demux.sif docker://devjune26/donor_demux:1.0

# r_eval container
singularity pull r_eval.sif docker://devjune26/r_eval:1.0
```

## Input Files

### Full Pipeline Mode
1. **Run Info CSV** (`runs_csv`): Metadata for donor BAM files. The example runs_csv can be found in https://github.com/june-zhang-bioinfo/donor-demux-nf/tree/main/data
   ```csv
   run_name,bam_file,ID,LB,PL,PU,SM,sex
   A1,/path/to/patient1.bam,P1_ID,P1_LIB,ILLUMINA,P1_UNIT,P1,Male
   A1,/path/to/patient2.bam,P2_ID,P2_LIB,ILLUMINA,P2_UNIT,P2,Female
   ```

2. **Reference Genome** (`reference`): FASTA file (e.g., GRCh38)

3. **10x Data Directory** (`tenx_dir`): Structure:
   ```
   run_name_1/
   ├── possorted_genome_bam.bam  
   ├── possorted_genome_bam.bam.bai
   ├── barcodes.tsv
   └── filtered_feature_bc_matrix.h5 (or .rds)
   ```

### Existing VCF Mode
Set `skip_vcf_prep = true` and provide:
- `vcf_path`: Path to pre-computed VCF file
- `donor_ids`: Comma-separated donor IDs
- `sexes`: Comma-separated sexes (Male/Female)

## Usage

### Full Pipeline (with VCF preparation)
```bash
module load nextflow
export JAVA_HOME="/apps/lib/java/jdk-20.0.2"
export PATH="$JAVA_HOME/bin:$PATH"
unset JAVA_CMD

nextflow run main.nf \
  -profile hpc \
  --runs_csv /data/run_info.csv \
  --reference /data/GRCh38.primary_assembly.genome.fa \
  --tenx_dir /data/10x/run_1 \
  --demux_outdir results/demuxalot
```

### Use pre-existing VCFs
```bash
nextflow run main.nf \
    -profile hpc \
    --skip_vcf_prep true \
    --reference "/path/GRCh38.primary_assembly.genome.fa" \
    --vcf_path "/path/data/test_dataset.vcf" \
    --tenx_dir "/path/data/10x/test_run"
```

## Parameters

| Parameter           | Required | Default                                   | Description                        |
| ------------------- | -------- | ----------------------------------------- | ---------------------------------- |
| `runs_csv`          | Yes*     | `/data/run_info.csv`                      | CSV with participant BAM metadata  |
| `reference`         | Yes*     | `/data/GRCh38.primary_assembly.genome.fa` | Reference genome FASTA             |
| `tenx_dir`          | Yes      | `/data/10x/run_name`                      | Directory containing 10x BAM files |
| `demux_outdir`      | No       | `results/demuxalot`                       | Output directory for results       |
| `skip_vcf_prep`     | No       | `false`                                   | Skip VCF preparation (demo mode)   |
| `vcf_path`          | Yes†    | `/data/test_dataset.vcf`                  | VCF file                           |
| `donor_ids`.        | Yes†    | (see config)                              | Donor IDs                          |
| `sexes`             | Yes†      (see config)                              | Donor sexes                        |

* Required unless skip_vcf_prep = true
† Required only when skip_vcf_prep = true

## Output Files

```
results/demuxalot/
└── [run_name]/
    ├── filtered.vcf
    ├── demuxalot_posteriors.csv       # Cell-to-donor assignment probabilities
    ├── demuxalot_likelihoods.csv      # Likelihood scores
    ├── PCA_plots.pdf        # PCA colored by sex and assignment
    └── SexGene_plots.pdf    # Sex-specific gene expression (XIST, RPS4Y1, etc.)
```


## Tips and troubleshooting

1. **Reuse vcf files**
   - If multiple experiments share the same donors, use the full pipeline first and get the `filtered.vcf`.
   - Use the `filtered.vcf` in vcf-skipping mode for other experiments.

2. **Clouds are not separated well in PC plots**
   - Check the number of SNPs in `filtered.vcf`.
  ```bash
  bcftools view -v snps filtered.vcf | wc -l
  ```
  - Consider to get better donor SNPs if it's below 100.

3. **"No RDS or H5 file found"**
   - Ensure 10x data directory contains either `.rds` or `filtered_feature_bc_matrix.h5`

4. **Memory issues**
   - Increase memory allocation in `nextflow.config`
   - Use `hpc` profile for larger resources

## Citation

If you use this pipeline, please cite:
- **Demuxalot**: Rogozhnikov, A., Ramkumar, P., Shah, K., Bedi, R., Kato, S., & Escola, G. S. *Demuxalot: scaled up genetic demultiplexing for single-cell sequencing*. bioRxiv (2021). https://doi.org/10.1101/2021.05.22.443646

- **FreeBayes**: Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012. https://arxiv.org/abs/1207.3907 

## Contact
	
Maintainer: June Zhang (zhang24925@gmail.com)
