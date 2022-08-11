#!/usr/bin/env python

# import numpy as np 
import bioinfo as bioinfo 
# import math as math 
import argparse
# import matplotlib as plt 
import gzip 
import itertools

def get_args():
    parser = argparse.ArgumentParser(description="specifies file inputs, index cutoff")
    parser.add_argument("-i", "--index", help="name of index list", \
        required=True, type=str)
    parser.add_argument("-f1", "--filename1", help=" name of input file", \
        required=True, type=str)
    parser.add_argument("-f2", "--filename2", help=" name of input file", \
        required=True, type=str)
    parser.add_argument("-f3", "--filename3", help=" name of input file", \
        required=True, type=str)
    parser.add_argument("-f4", "--filename4", help=" name of input file", \
        required=True, type=str)
    parser.add_argument("-ic", "--index_cutoff", help="Minimum quality score of an index", \
        required=True, type=int)
    return parser.parse_args()
#     # gives CLI options for file inputs and outputs


#shorten variable name
clinput = get_args()
dex= clinput.index
f1= clinput.filename1
f2= clinput.filename2
f3= clinput.filename3
f4= clinput.filename4
ic= clinput.index_cutoff
print("done getting CL input")
#  ./att.py -i ../TEST-input_FASTQ/index_list.txt -f1 ../TEST-input_FASTQ/test_r1.fq.gz -f2 ../TEST-input_FASTQ/test_r2.fq.gz -f3 ../TEST-input_FASTQ/test_r3.fq.gz -f4 ../TEST-input_FASTQ/test_r4.fq.gz -ic 30
# f1 = "../TEST-input_FASTQ/test_r1.fq.gz"
# f2 = "../TEST-input_FASTQ/test_r2.fq.gz"
# f3 = "../TEST-input_FASTQ/test_r3.fq.gz"
# f4 = "../TEST-input_FASTQ/test_r4.fq.gz"
# index_list = "../../../../../shared/2017_sequencing/indexes.txt" #switch to absolute path/ arg parse
index_dict = {} #create empty dictionary for list of indexes
indexes = set() #transfrom the dictionary to a set to make it
i=0 #line counter

with open(dex, "r") as ind: #populates the indexes set with the index list from the index list file
    ind.readline() #skips header line of colu,mns
    for line in ind:
        i +=1
        x = line.split("\t") 
        indexes.add(x[-1].strip("\n")) #adds the last column next to the new line character of each line (index) to the set
    # print(indexes)

# Reverse Compliment Function
bases = "ACGTN"
comps = "TGCAN"
tab = str.maketrans(bases,comps)
def rev_comp(seq: str) -> str:
    return seq.translate(tab)[::-1]

# print(rev_comp("ATCGNT"))

# Makes a set of the reverse compliments of the OG index list using rev_comp function
rev_comps = set()
rev_comp_dict = {}
for index in indexes:
    reverse=rev_comp(index)
    rev_comps.add(reverse)
    # print(reverse, index)
    rev_comp_dict[index]= reverse


#dictionary for unknown pairs
unknown_dict = {"unknown": 0}

#The below code creates a set of the observed hopped index permutations
hopped_dict = {}
perm_set = set()
x=list(itertools.permutations(indexes,2)) #uses itertools so make a permutation for index swapping in pairs/duos
for duo in x:
    x=duo[0]+"+"+duo[1]
    perm_set.add(x)
    hopped_dict[x] = 0
# print(hopped_dict)

matches = {} 
for index in indexes: #this loop opens files for matched input to be written to
    matches[index] = [open('output/matched/matched_'+index+"R1.fq", "w"), open('output/matched/matched_'+index+"R2.fq", "w")]
# print(matches)

matches_dict = dict.fromkeys(indexes, 0) # makes index set values into keys ina dictionary

#creating output files to write to for hopped and unmatched reads
unknown = (open("output/unknown/unknown_R1.fq", 'w'), open("output/unknown/unknown_R2.fq", "w"))
hopped = (open("output/hopped/hopped_R1.fq", "w"), open("output/hopped/hopped_R2.fq", "w"))
Leslie_count=0

with gzip.open(f1, "rt") as r1, gzip.open(f2, "rt") as r2, gzip.open(f3, "rt") as r3, gzip.open(f4, "rt") as r4: #this will open our 4 input file ONCE
    while True:
        R1_record = [r1.readline().strip(),r1.readline().strip(),r1.readline().strip(),r1.readline().strip()] #saves each line of a record in a list assigned to a single variable
        if R1_record[0] == '':
            break #this loop break stops the loop when it reaches the end of the file, files end with a blank string
        R2_record = [r2.readline().strip(),r2.readline().strip(),r2.readline().strip(),r2.readline().strip()]
        R3_record = [r3.readline().strip(),r3.readline().strip(),r3.readline().strip(),r3.readline().strip()]
        R4_record = [r4.readline().strip(),r4.readline().strip(),r4.readline().strip(),r4.readline().strip()]
        Leslie_count+=1
        # print(R1_record)
        # if R2_record[1] is not in indexes:
        if "N" in R2_record[1] or "N" in R3_record[1]: #writes low quality indexes to unmatched files
            unknown[0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
            unknown[1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
            unknown_dict["unknown"]+=1
        elif R2_record[1] not in indexes or R3_record[1] not in rev_comps:
            unknown[0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
            unknown[1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
            unknown_dict["unknown"]+=1
        elif bioinfo.qual_score(R2_record[3]) < ic or bioinfo.qual_score(R3_record[3]) < ic:
            unknown[0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
            unknown[1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
            unknown_dict["unknown"]+=1
        else:         # if R2_record[1] in indexes:
            if rev_comp_dict[R2_record[1]]==R3_record[1]: #check to see if indexes match
                #if bioinfo.qual_score(R2_record[3]) >= ic: #quality score check
                matches[R2_record[1]][0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
                matches[R2_record[1]][1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
                matches_dict[R2_record[1]]+=1
                # else:
                #     unknown[0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
                #     unknown[1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
                #     unknown_dict["unknown"]+=1
            else:   #if rev_comp_dict[R2_record[1]] != R3_record[1]:
                #if R3_record[1] in rev_comps:
                #    if bioinfo.qual_score(R3_record[3]) >= ic:
                hopped[0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
                hopped[1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
                hopped_dict[R2_record[1]+"+"+rev_comp(R3_record[1])]+=1
                    # else:
                    #     unknown[0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
                    #     unknown[1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
                    #     unknown_dict["unknown"]+=1
        # else:
        #     unknown[0].write(R1_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R1_record[1]+"\n"+R1_record[2]+"\n"+ R1_record[3]+"\n")   
        #     unknown[1].write(R2_record[0]+"_"+R2_record[1]+"+"+rev_comp(R3_record[1])+"\n"+ R2_record[1]+"\n"+R2_record[2]+"\n"+ R2_record[3]+"\n")
        #     unknown_dict["unknown"]+=1

# total_reads = len(r1.readlines())
total_reads = (sum(matches_dict.values())+sum(hopped_dict.values())+sum(unknown_dict.values()))
print("Total Records")
print(total_reads, "calc by sum", Leslie_count, "count by line", sep="\t")
print()
print("File Type"+"\t"+"Frequency"+"\t"+"Percentage of Total Reads")
print("Matched"+"\t"+str(sum(matches_dict.values()))+"\t"+str(sum(matches_dict.values())/total_reads))
print("Hopped"+"\t"+str(sum(hopped_dict.values()))+"\t"+str(sum(hopped_dict.values())/total_reads))
print("Unknown"+"\t"+str(sum(unknown_dict.values()))+"\t"+str(sum(unknown_dict.values())/total_reads))
print()
print("Matched Pair"+"\t"+"Frequency"+"\t"+"Percentage of Total Reads")

for pair in matches_dict:
    print(pair,matches_dict[pair],matches_dict[pair]/total_reads,sep="\t")
print()
print("Hopped Pair"+"\t"+"Frequency"+"\t"+"Percentage of Total Reads")   
for hop in hopped_dict:
    print(hop,hopped_dict[hop],hopped_dict[hop]/total_reads,sep="\t")
# print(unknown_dict["unknown"])

        

        




