#!/usr/bin/env bash

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demulti_pt2  ### Job Name
#SBATCH --output=slurm-%j-%x.out         ### File in which to store job output
#SBATCH --time=0-36:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=8       ### Number of cpus (cores) per task
#SBATCH --mail-user=abbiefayeolson@gmail.com
#SBATCH --mail-type=ALL

conda activate afo_py3

cd /home/afo/bgmp/bi624/demulti_pt2

/usr/bin/time -v ./PT2.py -R1 1294_S1_L008_R1_001.fastq.gz -I1 ./1294_S1_L008_R2_001.fastq.gz -I2 ./1294_S1_L008_R3_001.fastq.gz -R2 ./1294_S1_L008_R4_001.fastq.gz \
    -b ./indexes.txt -q 30
