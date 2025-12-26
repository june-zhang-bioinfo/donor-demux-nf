# Donor Demultiplexing Nextflow Pipeline

A Nextflow pipeline for demultiplexing pooled single-cell RNA-seq samples using genetic variants. The pipeline identifies individual donors in multiplexed 10x Genomics experiments through variant calling and probabilistic assignment.
This pipeline is intended for multiplexed 10x scRNA-seq experiments where donors are genetically distinct and no cell hashing (HTO) data is available.

## Overview

This pipeline performs:
1. **VCF Preparation** (optional): Processes donor BAM files to generate donor-specific variants
2. **Demultiplexing**: Uses Demuxalot to assign cells to donors based on genetic variants
3. **Quality Assessment**: Generates PCA plots and sex-specific gene expression visualizations

## Assumptions
- Donors are genetically unrelated (or weakly related)
- Adequate SNP coverage is present in scRNA-seq BAMs


## Pipeline Architecture

### Main Workflow (`main.nf`)
- Orchestrates the entire pipeline
- Supports two modes: full pipeline or demo mode with pre-existing VCF

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
- **Nextflow** (≥22.0)
- **Docker** or **Singularity**

### Container Images
- `devjune26/donor_demux:1.0` - Contains samtools, bcftools, FreeBayes, and demuxalot
- `devjune26/r_eval:1.0` - Contains R with Seurat and visualization packages

## Input Files

### Full Pipeline Mode
1. **Run Info CSV** (`params.runs_csv`): Metadata for donor BAM files
   ```csv
   run_name,bam_file,ID,LB,PL,PU,SM,sex
   A1,/path/to/patient1.bam,P1_ID,P1_LIB,ILLUMINA,P1_UNIT,P1,Male
   A1,/path/to/patient2.bam,P2_ID,P2_LIB,ILLUMINA,P2_UNIT,P2,Female
   ```

2. **Reference Genome** (`params.reference`): FASTA file (e.g., GRCh38)

3. **10x Data Directory** (`params.tenx_dir`): Structure:
   ```
   10x/
   ├── run_name_1/
   │   ├── possorted_genome_bam.bam
   │   ├── barcodes.tsv
   │   └── filtered_feature_bc_matrix.h5 (or .rds)
   └── run_name_2/
       └── ...
   ```

### Demo Mode
Set `params.skip_vcf_prep = true` and provide:
- `params.demo_vcf_path`: Path to pre-computed VCF file
- `params.demo_sample_names`: Comma-separated donor IDs
- `params.demo_sexes`: Comma-separated sexes (Male/Female)

## Usage

### Full Pipeline (with VCF preparation)
```bash
nextflow run main.nf \
  --runs_csv /data/run_info.csv \
  --reference /data/GRCh38.primary_assembly.genome.fa \
  --tenx_dir /data/10x \
  --demux_outdir results/demuxalot
```

### Demo Mode (skip VCF preparation)
```bash
nextflow run main.nf \
  --skip_vcf_prep \
  --demo_vcf_path /data/donors.vcf \
  --demo_sample_names "113_113,349_350,352_353" \
  --demo_sexes "Male,Male,Female" \
  --tenx_dir /data/10x \
  --demux_outdir results/demuxalot
```

### Using Profiles
```bash
# Local execution
nextflow run main.nf -profile local

# HPC/Slurm cluster
nextflow run main.nf -profile hpc
```

## Parameters

| Parameter           | Required | Default                                   | Description                        |
| ------------------- | -------- | ----------------------------------------- | ---------------------------------- |
| `runs_csv`          | Yes*     | `/data/run_info.csv`                      | CSV with participant BAM metadata  |
| `reference`         | Yes*     | `/data/GRCh38.primary_assembly.genome.fa` | Reference genome FASTA             |
| `tenx_dir`          | Yes      | `/data/10x`                               | Directory containing 10x BAM files |
| `demux_outdir`      | No       | `results/demuxalot`                       | Output directory for results       |
| `skip_vcf_prep`     | No       | `false`                                   | Skip VCF preparation (demo mode)   |
| `demo_vcf_path`     | Yes†     | `/data/text_dataset.vcf`                  | VCF file for demo mode             |
| `demo_sample_names` | Yes†     | (see config)                              | Donor IDs for demo mode            |
| `demo_sexes`        | Yes†     | (see config)                              | Donor sexes for demo mode          |

* Required unless skip_vcf_prep = true
† Required only when skip_vcf_prep = true

## Output Files

```
results/demuxalot/
└── [run_name]/
    ├── [run_name]_posteriors.csv       # Cell-to-donor assignment probabilities
    ├── [run_name]_likelihoods.csv      # Likelihood scores
    ├── [run_name]_PCA_plots.pdf        # PCA colored by sex and assignment
    └── [run_name]_SexGene_plots.pdf    # Sex-specific gene expression (XIST, RPS4Y1, etc.)
```

## Configuration Profiles

### `local`
- Single machine execution
- 2 CPUs, 6 GB memory per process
- 1 hour time limit

### `hpc`
- Slurm cluster execution
- 4 CPUs, 16 GB memory (default)
- FreeBayes: 8 CPUs, 32 GB, 12 hours
- R analysis: 1 CPU, 16 GB
- Automatic retry on failure (2 attempts)

## Building Docker Images

```bash
module load apptainer
# donor_demux container
apptainer pull donor_demux.sif docker://devjune26/donor_demux:1.0
# r_eval container
apptainer pull r_eval.sif docker://devjune26/r_eval:1.0
```

## Troubleshooting

### Common Issues

1. **"Run info CSV file not found"**
   - Ensure `params.runs_csv` points to a valid file
   - Or use `--skip_vcf_prep` for demo mode

2. **"BAM file not indexed"**
   - Pipeline automatically indexes BAM files if needed

3. **"No RDS or H5 file found"**
   - Ensure 10x data directory contains either `.rds` or `filtered_feature_bc_matrix.h5`

4. **Memory issues**
   - Increase memory allocation in `nextflow.config`
   - Use `hpc` profile for larger resources

## Citation

If you use this pipeline, please cite:
- **Demuxalot**: Heaton et al. (2020)
- **FreeBayes**: Garrison & Marth (2012)
- **Seurat**: Hao et al. (2021)

## Contact
	
Maintainer: June Zhang (zhang24925@gmail.com)
