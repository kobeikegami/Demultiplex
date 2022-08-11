
# Author: Kobe Ikegami email: kikegami@uoregon.edu
# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - 
#https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module
'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.6"         



def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return int(ord(letter)-33)


def qual_score(phred_score: str) -> float:
    """a loop to iterate over each letter in the phred_score and print them out.
    prints the letter and the corresponding score. Be sure to use your convert_phred function too."""
    i = 0
    sum = 0
    for letter in phred_score:
        s = convert_phred(letter)
        i += 1
        sum = sum + s
    return(sum / len(phred_score))  


# def validate_DNA_seq(DNA):
#     '''This function takes a string. Returns True if string is composed
#     of only As, Ts, Gs, and Cs. False otherwise. Case insensitive.'''
#     DNA = DNA.upper()
#     return len(DNA) == DNA.count("A") + DNA.count("T") + DNA.count("G") + DNA.count("C")

# assert validate_DNA_seq("aaaaa") == True, "DNA string not recognized"
# print("Correctly identified a DNA string")
# assert validate_DNA_seq("Hi there!") == False, "Non-DNA identified as DNA"
# print("Correctly determined non-DNA")
    
# DNA_bases = "ATCG"
# RNA_bases = "AUGC"
# def gc_content(DNA):
#     '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
#     DNA = DNA.upper()         #Make sure sequence is all uppercase
#     Gs = DNA.count("G")       #count the number of Gs
#     Cs = DNA.count("C")       #count the number of Cs
#     return (Gs+Cs)/len(DNA)

# assert gc_content("GCGCGC") == 1, "messed up calc when all GC"
# assert gc_content("AATTATA") == 0
# assert gc_content("GCATGCAT") == 0.5
# print("correctly calculated GC content")
    

# def oneline_fasta():
#     '''takes header line of FASTA file, isolates Gene name (key), 
#     protein ID and gene name (values) in a dicitonary. Then stores a protein sequence broken up
#     over multiple lines and concerts it into one line, saves it as a value'''
#  #!/usr/bin/env python
# import re
# # file"<filename>"

# gene_dict={}

# with open(homo_fa, "r") as fh:
#     protein='' #initialize empty variable
#     gene_id='' #""
#     for line in fh:
#         seq=line.strip('\n') #strip new line from header
#         if line.startswith(">"): #targets header by using the >
            
#             prot_id=re.findall(r"(?<=>)[a-z, A-Z, 0-9]+", line)[0] #RegEx to find prot ID and save as a variable
#             gene_id=re.findall(r"(?<=gene:)[a-z, A-Z, 0-9]+", line)[0] #RegEx to find gene ID and save as a variable
#             try:
#                 gene_name=re.findall(r"(?<=gene_symbol:)[a-z, A-Z, 0-9]+\s", line)[0] #RegEx to find gene name and save as variable
#             except:
#                 gene_name="" #makes exception not to grab gene name if there is none
        
#             protein="" # empty variable for next go round
#         else:
#             protein += line.strip("\n")
#         if gene_id in gene_dict:
#             if len(protein) > len(gene_dict[gene_id][0]):
#                 gene_dict[gene_id] = [protein, prot_id, gene_name]
#         else:
#             gene_dict[gene_id] = [protein, prot_id, gene_name]
# # print(len(gene_dict))

# def init_list(lst: list, value: float=0.0) -> list:
#     '''This function takes an empty list and will populate it with
#     the value passed in "value". If no value is passed, initializes list
#     with 101 values of 0.0.'''
#     lane1 = "../../lane1_NoIndex_L001_R1_003.fastq"
#     for x in range(0,101):
#         lst.append(value)
#     #print(lst)
#     return lst

# def populate_list(file: str) -> tuple[list, int]:
#     """creates a list of the sum of each qscore at each position of all sequences in a file """
    
#     test_list = [] #creates a list
#     test_list = init_list(test_list) #empty list to be populated with values passed as value, if not given a value, will initialize to 0
    
#     with open(file,"r") as fh:
#         i = 0 #counter for the entire file
#         for line in fh:
# #             print(i)
#             line = line.strip('\n')
#             i+=1
# #             print(line)
#             if i%4 == 0: #strips every 4th line (phred score line) 
#                 phred_score = line
# #                 print(phred_score)
#                 for count, ps in enumerate(phred_score): 
#                     test_list[count] += bioinfo.convert_phred(ps)
# #                     print(test_list)
#         return(test_list, i)

# bases = "ACGTN"
# comps = "TGCAN"
# tab = str.maketrans(bases,comps)
# def rev_comp(seq: str) -> str:
#     return seq.translate(tab)[::-1]

# print(rev_comp(""))