#!/bin/sh
#SBATCH --time 4-00:00:00
#SBATCH -p RT -n 700
module load mpi/mpich-3.2-x86_64

mpirun ./a.out
