#!/usr/bin/env python
import re

def softClipper(cigar:str, strand:str, position:int) -> int:
    ''' takes a cigar string, strand, and position and calculates the 5' start position'''
    if strand not in ("+","-"):
        print("Invalid entry: strand variable must be set to either + or -.")
        
    elif type(position) != int:
        print ("Invalid entry: position variable must be an integer.")
        
    elif re.fullmatch(r"([\d]+[MIDSN])+", cigar) == None:
        print ("Invalid entry: CIGAR string must have a number followed by chrs MIDSN to have a valid structure.")
    else:
        if strand == "+":
            leftS = re.findall(r"^(\d+)S", cigar)
            if len(leftS) > 0: # if there are any 'S' in cigar string:
                #print(leftS)
                leftS = leftS[0] #only take the leftmost S
                modify = int(leftS)
            else: # if there are no 'S' in cigar string
                modify = 0  
            return position - modify 
        else: # on negative strand
            # in any cigar string, there are at maximum 2 S's - one for each end of the string. 
            # first, sum up everything in the middle of the cigar string that we care about: M/D/N 
            mdn = re.findall(r"(\d+)[MDN]", cigar)
            #print(mdn)
            mdn_int = [int(i) for i in mdn]
            modify = sum(mdn_int)

            if re.search(r"(\d+)S$", cigar):
                rightS = re.findall(r"(\d+)S$", cigar)
                #print(rightS)
                modify += int(rightS[0])
                #print(modify)

            return position + modify 
        

#print(softClipper("haha", "+", 10))
#print(softClipper("hardeeharharhar", "-", 3923959))
#print(softClipper("10M38D283N", "-", 39234.58))
#print(softClipper("10M38D283N", "*", 12344))
#print(softClipper("10M38D283N", "+", 12344))
#print(softClipper("10M38D283N", "-", 12344))
#print(softClipper("2S8M", "+", 10))
#print(softClipper("2S8M", "-", 10))
#print(softClipper("2S8M2S", "-", 10))
# 20
#print(softClipper("1S27M113N38M", "-", 10))
# 188
#print(softClipper("36M12071N29M1S", "-", 10))
# 12147
#print(softClipper("58M8S", "-", 10))
# 76
#print(softClipper("8M2S", "-", 555))
# 565
#print(softClipper("10M", "-", 555))
# 565
#print(softClipper("2S8M", "+", 10))
#print(softClipper("8M5S", "+", 10))

#examples:
#softClipper(10M, 0, 10) == 10
#softClipper(2S7M, 0, 10) == 8
#softClipper(10M, 16, 10) == 19
#softClipper(10M5S, 16, 10) == 14

def strandParser(flag:int) -> str:
    ''' takes FLAG column 2 in sam file and interprets which strand the read is from '''
    if ((flag & 16) == 16):
        strand = "-"
    else:
        strand = "+"
    return strand

# example
#print(strandParser(0))
# +
#print(strandParser(16))
# -
#print(strandParser(256))
# +
#print(strandParser(291))
# +
#print(strandParser(307))
# -
#print(strandParser(1075))
# -

def qnameParser(qname:str) -> str:
    ''' takes qname string and parses it to return the umi '''
    return qname.split(":")[-1]

# example
#print(qnameParser("NS500451:154:HWKTMBGXX:1:11101:24260:1121:AAA"))
# AAA
#print(qnameParser("NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC"))
# CTGTTCAC
#print(qnameParser("11101:10568:1142:GAGAAGTC"))
# GAGAAGTC