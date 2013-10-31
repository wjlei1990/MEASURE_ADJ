#!/bin/bash
# parallel job using 
#PBS -l nodes=1:ppn=16,walltime=00:25:00

cd /tigress/lei/WORKFLOW7

num_nodes=6

module load openmpi
mpiexec -n $num_nodes ./TEST
