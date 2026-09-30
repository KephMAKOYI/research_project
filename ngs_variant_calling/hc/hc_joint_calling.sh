#!/bin/bash
#SBATCH --job-name=HC_ARDI
#SBATCH --account=b16-cbio-ag            # REPLACE with your project account
#SBATCH --partition=Main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1                # Head job needs minimal CPUs
#SBATCH --mem=8GB                        # Head job needs minimal RAM
#SBATCH --time=7-00:00:00                # Sufficient time for 10x 26GB WGS samples
#SBATCH --output=logs/nxf_head_%j.out
#SBATCH --error=logs/nxf_head_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kmakoyimakoyi@gmail.com       # REPLACE with your email

# 1. Load required environments
module load nextflow
module load singularity
module load samtools/1.19
module load bcftools/1.22
module load gatk/4.4.0.0
module load vep/106.1
module load htslib/1.19.1


# 2. Set directory paths (Using scratch3 for large intermediate WGS files)
# IMPORTANT: Replace <your_username> with your actual ilifu username
export NXF_SINGULARITY_CACHEDIR="/scratch3/users/$USER/singularity_cache"
export HC_WORK_DIR="/scratch3/users/$USER/work_dir/hc_joint_call/work_$(date +%Y%m%d)"

mkdir -p logs
mkdir -p $NXF_SINGULARITY_CACHEDIR

# 3. Configure Nextflow environment for Slurm stability
export NXF_OPTS="-Xms2G -Xmx6G"

# 4. Run Nextflow HC_JOINT_CALLING
# -c nextflow.config explicitly links to your custom configuration file
nextflow run /cbio/users/$USER/varcall-dsl2-2026/hc/hc_joint_calling.nf \
    -c /cbio/users/$USER/varcall-dsl2-2026/hc/hc_joint_calling.config \
    -work-dir $HC_WORK_DIR \
    --input /cbio/users/$USER/varcall-dsl2-2026/hc/hc_samplesheet.csv \
    -resume \
    -with-report logs/report_$(date +%F_%H-%M-%S).html \
    -with-timeline logs/timeline_$(date +%F_%H-%M-%S).html
