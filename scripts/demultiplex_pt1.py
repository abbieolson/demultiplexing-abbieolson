#!/usr/bin/env python

# PART 1

# read1, 1294_S1_L008_R1_001.fastq.gz
# index1, 1294_S1_L008_R2_001.fastq.gz
# index2, 1294_S1_L008_R3_001.fastq.gz
# read2, 1294_S1_L008_R4_001.fastq.gz

import argparse
import gzip
import numpy as np
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser(description="Program for demultiplexing")
parser.add_argument("-f", "--file", help="", required=True, type=str)

args = parser.parse_args()

FILE = args.file

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter) - 33

q_array = []
mean_scores = []
num_reads = 0
with gzip.open(FILE, "rt") as fh:
    LN = 0
    for line in fh:
        LN += 1
        line = line.strip()
        if LN == 4:
            for num in range(len(line)):
                q_array.append(0.0)
        if LN % 4 == 0:
            num_reads += 1
            for num in range(len(line)):
                #sys.exit()
                q_array[num] += convert_phred(line[num])

for item in q_array:
    mean_scores.append(item/num_reads)

plt.bar(range(0, len(mean_scores)), height = mean_scores, width=0.5)
plt.xlabel("mean quality score")
plt.ylabel("read position")
plt.title(os.path.basename(FILE))
plt.savefig(os.path.basename(FILE)+".png")


# GOOD Q-SCORE CUT-OFF: 30
