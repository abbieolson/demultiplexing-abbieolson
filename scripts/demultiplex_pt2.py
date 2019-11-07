#!/usr/bin/env python

# Ackowledgements:
# "While loop" and validation methodology inspired by Jared Galloway and Thomas Biondi
# File summary and dynamic out file inspiration from Wyatt Eng
# Reverse Complement inspiration from Stack Overflow
# Back-and-forth on keys and dictionaries with Pete Batzel
# Discussion on barcode incrementer with Anastasiya Prymolenna

import argparse
import gzip
import numpy as np

parser = argparse.ArgumentParser(description="Program for demultiplexing")

parser.add_argument("-R1", "--R1", help="", required=True, type=str)
parser.add_argument("-I1", "--I1", help="", required=True, type=str)
parser.add_argument("-I2", "--I2", help="", required=True, type=str)
parser.add_argument("-R2", "--R2", help="", required=True, type=str)
parser.add_argument("-b", "--barcode", help="", required=True, type=str)
parser.add_argument("-q", "--qscore", help="", required=True, type=int)

args = parser.parse_args()

R1 = args.R1
I1 = args.I1
I2 = args.I2
R2 = args.R2
barcode_file = args.barcode
qscore_cutoff = args.qscore

# set a global variable for reverse complementing
all_nucleotides = {'A':'T','C':'G','G':'C','T':'A','N':'N'}

def convert_phred(letter):
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def mean_qscore(phred):
    '''Converts string of quality scores to integers and finds the average'''
    score_array = []
    avg = 0

    for i in range(len(phred)):
        score_array.append(0)
    increment = 0
    for letter in phred:
        score_array[increment] += (convert_phred(letter))
        increment+=1
    mean_scores = np.array(score_array).mean()
    return mean_scores

def reverse_complement(seq):
    '''This function takes a sequence and returns its reverse complement'''
    for base in reversed(seq):
        return "".join(all_nucleotides[base] for base in reversed(seq))

def valid_read(read_index, barcodes, phred, q_cutoff):
    barcode_detected = read_index in barcodes
    read_passes = mean_qscore(phred) > q_cutoff
    return barcode_detected & read_passes

 # read in the index file and initialize a dictionary
barcodes_fi = barcode_file

raw_bcs = []
index_dict = {}
barcode_dict = {}

hopped_reads = 0
unknown_reads = 0
valid_reads = 0
total_reads = 0

# hopped reads
index_dict['R1_hopped'] = open('R1_hopped.fq', 'w')
index_dict['R2_hopped'] = open('R2_hopped.fq', 'w')

# unknown reads
index_dict['R1_unknown'] = open('R1_unknown.fq', 'w')
index_dict['R2_unknown'] = open('R2_unknown.fq', 'w')

barcodes_raw=open(barcodes_fi,"r")
# use readline() to remove header
barcodes_raw.readline()
lines = barcodes_raw.readlines()
for x in lines:
    barcode = x.strip().split()[4]
    raw_bcs.append(barcode)

    # dual matched
    index_dict['R1_' + barcode] = open('R1_' + barcode + '.fq', 'w')
    index_dict['R2_' + barcode] = open('R2_' + barcode + '.fq', 'w')

all_barcodes = raw_bcs

read_1 = gzip.open(R1, 'rt')
index_1 = gzip.open(I1, 'rt')
index_2 = gzip.open(I2, 'rt')
read_2 = gzip.open(R2, 'rt')

# 4 lists for each file for current record (length of 4)
# append 4 lines of record to list
while True:

    R1_record = []
    I1_record = []
    I2_record = []
    R2_record = []

    for x in range(4):
        R1_record.append(read_1.readline().strip())
        I1_record.append(index_1.readline().strip())
        I2_record.append(index_2.readline().strip())
        R2_record.append(read_2.readline().strip())

    # break out of loop if first line is an empty string
    if R1_record[0] == '':
        break
    if I1_record[0] == '':
        break
    if I2_record[0] == '':
        break
    if R2_record[0] == '':
        break

    # get indices and quality scores for each read
    # the index for read 2 must be reverse complemented
    read_1_index = I1_record[1]
    read_1_qscore = R1_record[3]
    read_2_index_rev = reverse_complement(I2_record[1])
    read_2_index = I2_record[1]
    read_2_qscore = R2_record[3]

    # check to see if the indices, barcodes, and quality scores return True
    valid_r1 = valid_read(read_1_index, all_barcodes, read_1_qscore, qscore_cutoff)
    valid_r2 = valid_read(read_2_index_rev, all_barcodes, read_2_qscore, qscore_cutoff)

    # update the headers in the output files to reflect the index pairs
    updated_r1_header = R1_record[0] + '_' + read_1_index + '_' + read_2_index
    updated_r2_header = R2_record[0] + '_' + read_1_index + '_' + read_2_index

    # if the input files pass the valid_read() function, pass into loop
    # writes 4-lined fastq files and modifies the headers

    read_concat = read_1_index + '_' + read_2_index

    if valid_r1 & valid_r2:
        # for dual-matched reads
        if read_1_index == read_2_index_rev:

            if read_concat not in barcode_dict:
                barcode_dict[read_concat] = 1
            else:
                barcode_dict[read_concat] += 1

            index_dict['R1_' + read_1_index].write(updated_r1_header)
            index_dict['R1_' + read_1_index].write(R1_record[0] + '_' + read_1_index + '_' + read_2_index + '\n')
            index_dict['R1_' + read_1_index].write(R1_record[1] + '\n')
            index_dict['R1_' + read_1_index].write(R1_record[2] + '\n')
            index_dict['R1_' + read_1_index].write(R1_record[3] + '\n')

            index_dict['R2_' + read_2_index_rev].write(updated_r2_header)
            index_dict['R2_' + read_2_index_rev].write(R2_record[0] + '_' + read_1_index + '_' + read_2_index + '\n')
            index_dict['R2_' + read_2_index_rev].write(R2_record[1] + '\n')
            index_dict['R2_' + read_2_index_rev].write(R2_record[2] + '\n')
            index_dict['R2_' + read_2_index_rev].write(R2_record[3] + '\n')
            valid_reads += 1
            total_reads += 1

        else:
            # for hopped reads
            if 'hopped' not in barcode_dict:
                barcode_dict['Hopped Reads'] = 1
            else:
                barcode_dict['Hopped Reads'] += 1

            index_dict['R1_hopped'].write(updated_r1_header)
            index_dict['R1_hopped'].write(R1_record[0] + '_' + read_1_index + '_' + read_2_index + '\n')
            index_dict['R1_hopped'].write(R1_record[1] + '\n')
            index_dict['R1_hopped'].write(R1_record[2] + '\n')
            index_dict['R1_hopped'].write(R1_record[3] + '\n')

            index_dict['R2_hopped'].write(updated_r2_header)
            index_dict['R2_hopped'].write(R2_record[0] + '_' + read_1_index + '_' + read_2_index + '\n')
            index_dict['R2_hopped'].write(R2_record[1] + '\n')
            index_dict['R2_hopped'].write(R2_record[2] + '\n')
            index_dict['R2_hopped'].write(R2_record[3] + '\n')
            hopped_reads += 1
            total_reads += 1

    else:
        # for unknown reads
        if 'unknown' not in barcode_dict:
            barcode_dict['Unknown Reads'] = 1
        else:
            barcode_dict['Unknown Reads'] += 1
        index_dict['R1_unknown'].write(updated_r1_header)
        index_dict['R1_unknown'].write(R1_record[0] + '_' + read_1_index + '_' + read_2_index + '\n')
        index_dict['R1_unknown'].write(R1_record[1] + '\n')
        index_dict['R1_unknown'].write(R1_record[2] + '\n')
        index_dict['R1_unknown'].write(R1_record[3] + '\n')

        index_dict['R2_unknown'].write(updated_r2_header)
        index_dict['R2_unknown'].write(R2_record[0] + '_' + read_1_index + '_' + read_2_index + '\n')
        index_dict['R2_unknown'].write(R2_record[1] + '\n')
        index_dict['R2_unknown'].write(R2_record[2] + '\n')
        unknown_reads += 1
        total_reads += 1

# close file handles at end of file
read_1.close()
index_1.close()
index_2.close()
read_2.close()

# outputs for the summary report
print("Begin of Summary Report:")
print("Total Reads = {}".format(total_reads))
print("Total Valid Reads = {}".format(valid_reads))
print("Total Hopped Reads = {}".format(hopped_reads))
print("Total Unknown Reads = {}".format(unknown_reads))

# prints percentages of barcodes that occurred and the overall amount of index swapping that occurred
# won't work in py2, must use py3

for key in barcode_dict:
    print("The percentage of {} that occurred is {}.".format(key, barcode_dict[key] / valid_reads * 100), "%")

print("The overall amount of index swapping that occurred is {}".format((hopped_reads/total_reads) * 100), "%")

print("End of Summary Report")
