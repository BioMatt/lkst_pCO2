#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=def-dogfish
#SBATCH --mem=20000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-66

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load salmon/1.7.0

list=/home/biomatt/scratch/pCO2/trimmed_reads.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}"
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2}')
set -- $var

forward=$1
reverse=$2

echo ${forward}
echo ${reverse}

sample=$(basename ${forward} .R1.fq.gz | cut -d "." -f3)

echo ${sample}

salmon quant -i /home/biomatt/scratch/pCO2/salmon_index -l IU -1 ${forward} -2 ${reverse} --validateMappings --seqBias --gcBias -o /home/biomatt/scratch/pCO2/salmon_quant/${sample}
