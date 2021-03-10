#!/bin/bash
###########__________Run script__________#############
################ Hpc machine ################

#$ -cwd -V
#$ -l h_rt=48:00:00
#$ -t 1-2
export HPC_MODE=TRUE
TASK_ID=1
python3 main.py $TASK_ID
