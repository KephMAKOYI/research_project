#!/bin/bash
#SBATCH --job-name=6JFK_MFN2_mT_job
#SBATCH --time=5-00:00:00            # Wall-clock time (e.g., 5 days)
#SBATCH --mem=500G                   # Memory per node (e.g., 500 GB)
#SBATCH --cpus-per-task=32           # Number of CPUs per task
#SBATCH --output=%x_%j.out           # Standard output file
#SBATCH --error=%x_%j.err            # Standard error file
#SBATCH --partition=HighMem          # Or other partition like Devel, GPU, HighMem
#SBATCH --mail-user=kmakoyimakoyi@gmail.com
#SBATCH --mail-type=ALL

module load gromacs/2024.2

INPUTS="/cbio/users/kmakoyi/outdir/vcfs/joint/SnpEff/gromacs/inputs"

# 1. strip out the crystal waters (To delete the water molecules (residue "HOH" in the PDB file))
#grep -v HOH 6JFK_MFN2_mT_model_0.pdb > 6JFK_MFN2_mT_clean.pdb

# 2. Generate topology
#gmx pdb2gmx -f 6JFK_MFN2_mT_clean.pdb -o 6JFK_MFN2_mT_processed.gro -water tip3p
# type 6 press Enter (AMBER99SB-ILDN)

# 3.define the box
#gmx editconf -f 6JFK_MFN2_mT_processed.gro -o 6JFK_MFN2_mT_newbox.gro -c -d 1.0 -bt cubic

# 4. fill box with solvent (water)
#gmx solvate -cp 6JFK_MFN2_mT_newbox.gro -cs spc216.gro -o 6JFK_MFN2_mT_solv.gro -p topol.top

# 5. Generate ions.mdp file
#gmx grompp \
 #   -f $INPUTS/ions.mdp \
 #   -c 6JFK_MFN2_mT_solv.gro \
 #   -p topol.top \
 #   -o 6JFK_MFN2_mT_ions.tpr

#gmx genion \
 #   -s 6JFK_MFN2_mT_ions.tpr \
 #   -o 6JFK_MFN2_mT_solv_ions.gro \
 #   -p topol.top -pname NA -nname CL -neutral
#choose group 13 "SOL"

# 6. Energy Minimization

# Remove steric clashes
# Generate minim.mdp

#gmx grompp \
 #   -f $INPUTS/minim.mdp \
 #   -c 6JFK_MFN2_mT_solv_ions.gro \
 #   -p topol.top \
 #   -o 6JFK_MFN2_mT_em.tpr

# Run
#gmx mdrun -v -deffnm 6JFK_MFN2_mT_em

# 7. Equilibration is often conducted in two phases.
# The first phase is conducted under an NVT ensemble (constant Number of particles, Volume, and Temperature).
# This ensemble is also referred to as "isothermal-isochoric" or "canonical."

# Get the nvt.mdp file
#gmx grompp \
 #   -f $INPUTS/nvt.mdp \
 #   -c 6JFK_MFN2_mT_em.gro \
 #   -r 6JFK_MFN2_mT_em.gro \
 #   -p topol.top \
 #   -o 6JFK_MFN2_mT_nvt.tpr

#gmx mdrun -deffnm 6JFK_MFN2_mT_nvt -v


# 8. NVT equilibration
#gmx grompp \
 #   -f $INPUTS/npt.mdp \
 #   -c 6JFK_MFN2_mT_nvt.gro \
 #   -r 6JFK_MFN2_mT_nvt.gro \
 #   -t 6JFK_MFN2_mT_nvt.cpt \
 #   -p topol.top \
 #   -o 6JFK_MFN2_mT_npt.tpr

#gmx mdrun -deffnm 6JFK_MFN2_mT_npt -v


# 9. Production MD
# Generate md.mdp
# We will run 0.1ns or 100ps

gmx grompp \
    -f $INPUTS/md.mdp \
    -c 6JFK_MFN2_mT_npt.gro \
    -t 6JFK_MFN2_mT_npt.cpt \
    -p topol.top \
    -o 6JFK_MFN2_mT_md.tpr

gmx mdrun -deffnm 6JFK_MFN2_mT_md -v

# Running GROMACS on GPU
#gmx mdrun -deffnm md_0_1 -nb gpu

