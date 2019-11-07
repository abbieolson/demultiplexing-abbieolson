#!/usr/bin/env bash

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demulti_pt2  ### Job Name
#SBATCH --output=slurm-%j-%x.out         ### File in which to store job output
#SBATCH --time=0-24:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=8       ### Number of cpus (cores) per task
#SBATCH --mail-user=abbiefayeolson@gmail.com
#SBATCH --mail-type=ALL

/usr/bin/time -v /home/afo/bgmp/bi624/demulti_pt2/PT2.py -R1 /home/afo/bgmp/bi624/demulti_pt2/demulti_pt1/demultiplexing-abbieolson/r1test.fastq.gz -I1 /home/afo/bgmp/bi624/demulti_pt2/demulti_pt1/demultiplexing-abbieolson/r2test.fastq.gz -I2 /home/afo/bgmp/bi624/demulti_pt2/demulti_pt1/demultiplexing-abbieolson/r3test.fastq.gz -R2 /home/afo/bgmp/bi624/demulti_pt2/demulti_pt1/demultiplexing-abbieolson/r4test.fastq.gz \
    -b /home/afo/bgmp/bi624/demulti_pt2/indexes.txt -q 30

