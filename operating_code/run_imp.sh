#!/bin/bash

#$ -q shai.q@bhn27
#$ -cwd
#$ -N impute
#$ -j yes

LD_LIBRARY_PATH=/gpfs0/shai/projects/software/R4/R-4.4.1/lib64/R/lib/:$LD_LIBRARY_PATH
/gpfs0/shai/projects/software/R4/R-4.4.1/bin/Rscript softimpute_bootstrapped.R $1 $2