#!/bin/bash
#MSUB -A b1042
#MSUB -l walltime=4:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -j oe
#MSUB -q genomics



module load shapeit
module load plink
module load vcftools
module load python
module load R



