#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=def-dogfish
#SBATCH --mem=20000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load salmon/1.7.0

salmon index -t /home/biomatt/scratch/pCO2/Trinity_gill.fasta -i /home/biomatt/scratch/pCO2/salmon_index -k 31