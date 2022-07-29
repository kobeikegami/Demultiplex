#!/bin/bash

#SBATCH --partition=bgmp       ### Partition (like a queue in PBS)
#SBATCH --job-name=atf_p1_r2    ### Job Name/optional
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node, optional
#SBATCH --account=bgmp          ### Account used for job submission/
#SBATCH --cpus-per-task=8

python atf_p1.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -l 101 -o "r2"