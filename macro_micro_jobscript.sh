#!/bin/bash
#
#SBATCH --job-name=macro-micro-kschidha
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#
#SBATCH --ntasks=5
#SBATCH --ntasks-per-node=4
#SBATCH --time=05:00:00
#
# Working directory:
#SBATCH -D ./

# load modules
module use /usr/local.nfs/sgs/modulefiles
module load vtk/9.0.1 cmake/3.18.2 ub2004/boost/1.75.0 ub2004/libxml2/2.9.10
# should also load gcc/10.2 openmpi/3.1.6-gcc-10.2 
# enable custom build openmpi 3.1 that works with slurm
export CPATH=/scratch-nfs/kschidha/openmpi/install-3.1/include
export PATH=/scratch-nfs/kschidha/openmpi/install-3.1/bin:$PATH

echo "Finished Loading"

# possibly add precice path

echo "Launching macro participant"
srun -n 1 macro-heat/test_macro_heat & 

echo "Launching micro manager"
srun -n 4 python3 python-micro-heat/run-micro-problems.py

echo "Simulation completed."



