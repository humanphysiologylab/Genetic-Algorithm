#!/bin/sh
#SBATCH --time 5:00:00
#SBATCH -p RT -n 8
#SBATCH -J TASK_NAME
#SBATCH --comment Atrial_genetic_algorithms

mpirun ./ga input_ga.txt &> out.log
