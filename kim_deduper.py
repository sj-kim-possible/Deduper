#!/usr/bin/env python
#
#.-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#/ / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
#`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'
# Deduper
#
# Bi624 2022
# sj kim
#
# script overview: Takes a sorted (via samtools) UNIQUELY MAPPED sam file and deduplicates - removes
# pcr duplicates and retains the first encountered duplicate. Outputs a sorted sam file.
# sorting command:
# samtools sort <input samfile> -o <output samfile>

#### [ port ] ####
import argparse
import math 
import re

def get_args():
    parser = argparse.ArgumentParser(description="This script removes PCR duplicates in a SORTED sam file of uniquely mapped reads.")
    parser.add_argument("-f", "--file", help="Sorted uniquely mapped input sam file.", required=True, type = str)
    parser.add_argument("-o", "--outfile", help="Sorted uniquely mapped DEDUPED output sam file.", required = True, type = str)
    parser.add_argument("-u", "--umi", help="txt file of Valid UMIs.", required=True, type = str)
    return parser.parse_args()
args = get_args()

#############################################################
#                                                           #
#                         deDuper                           #
#                                                           #
#############################################################

# import umis into a set
with open(args.umi) as umiFile:
    validUmis = {umi.strip('\n') for umi in umiFile.readlines()}

#print(validUmis)

##################################
#### [ high-level functions ] ####
##################################

def qnameParser(qname:str) -> str:
    ''' takes QNAME string and parses it to return the umi '''
    return qname.split(":")[-1]

def strandParser(flag:int) -> str:
    ''' takes FLAG column 2 in sam file and interprets which strand the read is from '''
    if ((flag & 16) == 16):
        strand = "-"
    else:
        strand = "+"
    return strand

def softClipper(cigar:str, strand:str, position:int) -> int:
    ''' takes CIGAR string, strand, and position and calculates the 5' start position'''
    if strand not in ("+","-"):
        print("Invalid entry: strand variable must be set to either + or -.")
        
    elif type(position) != int:
        print ("Invalid entry: position variable must be an integer.")
        
    # r"" represents that the string in re will be interpreted literally - escape characters won't be escaped
    elif re.fullmatch(r"([\d]+[MIDSN])+", cigar) == None:
        print ("Invalid entry: CIGAR string must have a number followed by chrs MIDSN to have a valid structure.")
    else:
        if strand == "+":
            leftS = re.findall(r"^(\d+)S", cigar)
            if len(leftS) > 0: # if there are any 'S' in cigar string:
                leftS = leftS[0] #grab the integer of leftS
                modify = int(leftS)
            else: # if there are no 'S' in cigar string
                modify = 0  
            return position - modify 
        else: # on negative strand
            # in any cigar string, there are at maximum 2 S's - one for each end of the string. 
            # first, sum up everything in the middle of the cigar string that we care about: M/D/N 
            mdn = re.findall(r"(\d+)[MDN]", cigar)
            mdn_int = [int(i) for i in mdn]
            modify = sum(mdn_int)

            if re.search(r"(\d+)S$", cigar):
                rightS = re.findall(r"(\d+)S$", cigar)
                modify += int(rightS[0])

            return position + modify 


#######################
#### [ de-duping ] ####
#######################

currChrom = ""                          #keeps track of current chromosome
uniqueRecordsByChrom = set()            #keeps track of unique records in tuple(pos, UMI, strand)

# counters
totalRecords = 0
uniqueRecords = 0
dupes = 0
invalidUmis = 0
headers = 0

# indexes
qname_index = 0
flag_index = 1
rname_index = 2
pos_index = 3
cigar_index = 5

with open(args.file, "r") as rawSamInfile, open(args.outfile, "w") as dedupedOutfile:
    while True:
        line = rawSamInfile.readline()
        if line == "":
            break
        #print(line)

# send header lines out
        if line.startswith("@"):
            dedupedOutfile.write(line)
            headers += 1
            continue

# turn each samfile line into a record list, separated by tabs
        record = line.strip("\n").split("\t") # record is a list
        totalRecords += 1

# if the record does not have a valid umi, continue to next record:
        umi = qnameParser(record[qname_index])
        if umi not in validUmis:
            invalidUmis += 1
            continue

# assuming valid umi, check strand and position:
        strand = strandParser(int(record[flag_index]))
        cigar = record[cigar_index]
        sam_position = int(record[pos_index])

        fivePrimeStartPos = softClipper(cigar, strand, sam_position)

# generate a tuple containing the identifiers for this record:
        identifierTup = (fivePrimeStartPos, umi, strand)
        chrom = record[rname_index]
        
# re-set the tuple set and chrom var for each new chromosome:
        if chrom != currChrom:
            uniqueRecordsByChrom.clear()
            currChrom = chrom

        if identifierTup not in uniqueRecordsByChrom:
            uniqueRecordsByChrom.add(identifierTup)
            uniqueRecords += 1
            dedupedOutfile.write(line)

        else: # identifierTup IS ALREADY IN unique records set, so it's a duplicate
            dupes += 1


# summary stats:
print(f"header lines: {headers}")
print(f"total records: {totalRecords}")
print(f"unique records: {uniqueRecords}")
print(f"duplicates: {dupes}")
print(f"invalid umis: {invalidUmis}")

# ./kim_deduper.py -f ./test.sam -u ./STL96.txt -o ./STL96_test.sam
# header lines: 24
# total records: 76
# unique records: 54
# duplicates: 21
# invalid umis: 1

# cat C1_SE_uniqAligned_deduped.sorted.sam | grep -v "^@" | awk '{print $3}' | uniq -c | sort -k2 -V
