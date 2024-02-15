#!/bin/bash
#SBATCH --job-name=cx
#SBATCH --partition=s_parallel
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=36
#SBATCH --time=2-00:00:00
#SBATCH --exclusive
module purge
module load tools/octopus/13.0

export OMP_NUM_THREADS=1
/cluster/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/openmpi-3.1.2-af7zrcesig4lnjlho2aewuwpozudwv2e/bin/orterun -np $SLURM_NTASKS octopus 
