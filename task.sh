#!/bin/sh
#SBATCH --time 2-0:00:00
#SBATCH -p RT -N 32 -n 512
#SBATCH -J ga_maleckar
#SBATCH --comment Atrial_genetic_algorithms

#module load mpi

mpirun ./ga input_ga.txt &> out.log
