#!/bin/bash
#SBATCH --job-name=BAM_NGS
#SBATCH --account=b16-cbio-ag          # REPLACE with your project account
#SBATCH --partition=HighMem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1              # Head job needs minimal CPUs
#SBATCH --mem=8GB                      # Head job needs minimal RAM
#SBATCH --time=11-00:00:00                # Sufficient time for 10x 26GB WGS samples
#SBATCH --output=logs/nxf_head_%j.out
#SBATCH --error=logs/nxf_head_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kmakoyimakoyi@gmail.com       # REPLACE with your email

# 1. Load required environments
module load nextflow
module load singularity
module load trimgalore/0.6.10
module load bwa/github
module load samtools/1.19
module load bcftools/1.22
module load gatk/4.4.0.0
module load htslib/1.19.1
module load vg/1.65.0

# 2. Set directory paths (Using scratch3 for large intermediate WGS files)
# IMPORTANT: Replace <your_username> with your actual ilifu username
export NXF_SINGULARITY_CACHEDIR="/scratch3/users/$USER/singularity_cache"
export WORK_DIR="/scratch3/users/$USER/work_dir/main/bam/work_$(date +%Y%m%d)"
export results="/cbio/users/$USER/results"

mkdir -p logs
mkdir -p $NXF_SINGULARITY_CACHEDIR
mkdir -p $results

# 3. Configure Nextflow environment for Slurm stability
export NXF_OPTS="-Xms2G -Xmx6G"

# 4. Run Nextflow
# -c nextflow.config explicitly links to your custom configuration file
nextflow run /cbio/users/$USER/varcall-dsl2-2026/main/bam/main_bam.nf \
    -c /cbio/users/$USER/varcall-dsl2-2026/main/bam/bam.config \
    -work-dir $WORK_DIR \
    --input /cbio/users/$USER/varcall-dsl2-2026/main/bam/bam_samplesheet.csv \
    --outdir $results \
    -resume \
    -with-report logs/report_$(date +%F_%H-%M-%S).html \
    -with-timeline logs/timeline_$(date +%F_%H-%M-%S).html
