#!/bin/bash
#SBATCH -p compute      # Select a partition to run.
#SBATCH -N 4           # Number of nodes
#SBATCH -t 06:00:00     # Request 40 hours walltime
#SBATCH -U HPC_12_00409 # Required project.

srun hostname

echo $PWD

export PSPOT_DIR=/home/users/mehese/CASTEP-6.0/pseudopotentials

mpirun -npernode 12 ~/CASTEP-6.0/obj/linux_x86_64_ifort11/castep tmplate
