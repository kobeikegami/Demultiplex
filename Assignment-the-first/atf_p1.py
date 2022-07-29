#!/usr/bin/env python

#script to make histograms of per-base distributiond of qualituy scores

import bioinfo
import gzip
import argparse
import matplotlib.pyplot as plt


#defines arguments to be input in CLI for script runs
def get_args():
    parser = argparse.ArgumentParser(description="specifies input, output files, and read length")
    parser.add_argument("-o", "--output", help=" specifies output file name", \
        required=True, type=str)
    parser.add_argument("-f", "--filename", help="name of input file", \
        required=True, type=str)
    parser.add_argument("-l", "--length", help="length of read", \
        required=True, type=int)
    return parser.parse_args()

#shorten variable names
clinput = get_args()
o= clinput.output
f= clinput.filename
l= clinput.length 

#absolute paths for all relevant fq files
# r1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
# i1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
# i2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
# r2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"     

num_lines=0  
def init_list(lst: list, value: float=0.0) -> list:
    for x in range(0,l):
        lst.append(value)
    #print(lst)
    return lst

test_list = [] #creates a list
test_list = init_list(test_list)
    
with gzip.open(f,"rt") as fh:
    i = 0 #counter for the entire file
    for line in fh:
        num_lines+=1
#             print(i)
        line = line.strip('\n')
        i+=1
#             print(line)
        if i%4 == 0: #strips every 4th line (phred score line) 
            phred_score = line
#                 print(phred_score)
            for count, ps in enumerate(phred_score): 
                test_list[count] += bioinfo.convert_phred(ps)
    # print(test_list)
    # return(test_list, i)

plot_list = init_list([])
for position, sums in enumerate(test_list):
    mean = (sums/((num_lines)/4))
    test_list[position] = mean
    plot_list[position] = position
# print(test_list)
# print(plot_list)

# def hist_plot(test_list, o, l):
#     x_list=[]
#     i=0
#     for i in range(0,l):
#         x_list.append(i)
#         i+=l

plt.bar(plot_list, test_list)
plt.xlabel("Base Position")
plt.ylabel("Mean Quality")
plt.title("Mean Quality Score Per Base position")
plt.savefig(o)
# plt.savefig("{}.png".format(o))