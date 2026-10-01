## Identification of novel rare disease genes and variants in the Democratic Republic of Congo using Machine Learning and in silico Protein Structure Prediction Approaches

Integrative bioinformatics pipeline combining GRCh38 and African Pangenome reference alignments with AI-driven pathogenicity prediction and molecular dynamics (MD) simulations to reclassify variants of uncertain significance (VUS) in rare diseases within the Democratic Republic of the Congo.

### Overview

This repository contains an end-to-end Nextflow workflow and MD simulation suite designed to improve diagnostic yields for African populations by:
1. Comparing traditional reference genome alignments against pangenomic graphs.
2. Integrating ensemble machine learning pathogenicity tools (AlphaMissense, CADD, REVEL, etc.).
3. Utilizing biophysical molecular dynamics simulations to evaluate structural impacts on candidate proteins.

![Pipeline Workflow](docs/assets/pipeline_workflow.PNG)

### Research Objectives
* Improved variant detection: Integrate the African pangenome graph into standard NGS workflows to enhance variant call precision in underrepresented populations.
* Improved variant Reclassification: Leverage machine learning ensemble scores and atomistic molecular dynamics to resolve Variants of Uncertain Significance (VUS).

### Methodology

1. Nextflow NGS Variant Calling
* Alignment: Standard mapping to GRCh38 via BWA-MEM alongside graph-based alignment to the African Pangenome via vg Giraffe.
* Variant Calling: Parallel benchmarking using GATK HaplotypeCaller and DeepVariant.
* Annotation & Prioritization: Functional annotation via SnpEff and dbNSFP scores (AlphaMissense, CADD, PolyPhen, REVEL, SIFT).

2. Structural & Molecular Dynamics (MD) Analysis
* System Preparation: Interactive structure building using CHARMM-GUI and AI-predicted 3D models.
* Simulation Engine: GROMACS execution evaluating wild-type (WT) vs. mutant structural stability across thermodynamic ensembles (NVT, NPT, Production MD).

### Repository Structure

![Pipeline Workflow](docs/assets/repository_structure.PNG)

### Getting Started

Prerequisites
* HPC Cluster: Slurm-managed HPC (tested on Ilifu HPC HighMem partition)
* Workflow Manager: Nextflow (>=22.04)
* Container / Software Engines: Docker / Singularity / Conda
* Simulation Suite: GROMACS (>=2021.x)

### Execution Example

To run the full NGS variant calling pipeline on a Slurm cluster:

cd ngs_variant_calling/main
sbatch nextflow.sh

To run a molecular dynamics simulation for a mutant protein structure:

cd molecular_dynamic_simulation
bash 6JFK_MFN2_mT.sh

### Ethical Considerations

* Confidentiality & Security: All human genomic data were processed under strict security protocols on the Ilifu HPC platform.
* Ethical Approval: Protocol approved by the Human Research Ethics Committee (HREC Ref No: 881/2025).

### Future Roadmap

* Expansion to larger and phenotypically diverse cohort sizes.
* In vitro experimental validation of top prioritized candidate variants.
* Systematic integration of structural variant (SV) discovery pipelines.
* Scaling MD trajectory production runs to $\ge 100\text{ ns}$ across all identified VUS candidates.

#### Authors & Affiliations

Authors: Keph Makoyi, Nicola Mulder, Christian D. Bope, Hocine Bendou, Aimé Lumaka
* University of Cape Town
* University of Kinshasa
* African Rare Disease Initiative
