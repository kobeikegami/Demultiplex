#!/usr/bin/bash
#SBATCH --partition=bgmp       ### Partition (like a queue in PBS)
#SBATCH --job-name=demux     ### Job Name/optional
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node, optional
#SBATCH --account=bgmp          ### Account used for job submission/
#SBATCH --cpus-per-task=1
#SBATCH --error=demux_%j.err
#SBATCH --out=demux_%j.out

conda activate base
mydir="/projects/bgmp/kikegami/bioinfo/Bi622/Demultiplex/Assignment-the-third/output"
datadir="/projects/bgmp/shared/2017_sequencing"

/usr/bin/time -v ./att.py -i $datadir/indexes.txt \
-f1 $datadir/1294_S1_L008_R1_001.fastq.gz \
-f2 $datadir/1294_S1_L008_R2_001.fastq.gz \
-f3 $datadir/1294_S1_L008_R3_001.fastq.gz \
-f4 $datadir/1294_S1_L008_R4_001.fastq.gz \
-ic 30


