#!/bin/bash
#
#SBATCH --job-name=VP85_430
#SBATCH --output=v4_1_85_angle_30.txt
#
#SBATCH -n 8
#SBATCH -p short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gibron@seznam.cz
mpirun python dip2.4.py