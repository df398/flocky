#!/bin/bash

#SBATCH --job-name=C60flocky
#SBATCH -A all-account
#SBATCH -p allq
#SBATCH --ntasks=16
#SBATCH --time=60:00:00
#SBATCH -e slurm.info.err
#SBATCH -o slurm.info.out

#On our system, we load the relevant compilers (typically those that were used to compile flocky, tapreaxff and associated tools)
source /opt/intelvars/compilervars_2019.sh intel64

#mpirun path may be modified according to local MPI environment. We also assume that flocky_mpi is set in your PATH.
mpirun -np $SLURM_NTASKS flocky_mpi

# End of submit file
