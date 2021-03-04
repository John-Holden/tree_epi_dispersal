#!/bin/bash
###########__________Run script__________#############
################ Hpc machine ################


SGE_TASK_ID=1
export HPC_MODE=True
python3 main.py SGE_TASK_ID
