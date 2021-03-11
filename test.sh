#!/bin/bash
###########__________Run script__________#############
################ Hpc machine ################

export HPC_MODE=TRUE
SGE_TASK_ID=1
python3 main.py $SGE_TASK_ID
