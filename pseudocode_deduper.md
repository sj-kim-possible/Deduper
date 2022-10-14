# Deduper

## Part 1: Pseudocode
Define the problem:
The overall goal of this script is to remove PCR-duplicate reads from a sam file using a reference genome. A PCR duplicate is an identical molecule made by PCR. PCR is used in sequencing to amplify species of interest so that there is enough DNA to sequence and also to quench the signal to noise ratio in sequencing. PCR duplicates create bias in the data by representing a sequence that isn't biologically present. Additionally, there is bias in amplification. For example, sequences with high G/C content doesn't amplify as well due to higher binding strength. Longer DNA sequences also don't amplify as easily. A sequence is a PCR duplicate when the chromosome, 5' start position, UMI, and strand are the same. Thus, this informs the strategy for removing PCR duplicates. The code must account for single-end data at minimum. Paired-end if I'm feelin it. Single-end data will result in more reads being flagged as duplicates since the data is less specific. 

argparse options for:
file: absolute path to sorted sam file that needs de-duplicating
output file: filepath for deduplicated sam file
umi: filepath for umi
help: prints a useful help message

*** high level functions ***
def softClipper(cigar (str), strand (str), position(int)): --> int
''' takes a cigar string, strand, and position and calculates the 5' start position'''
    if strand == +:
        fivePrimeStartPos = posStrandCigarParser(cigar, position)
    else strand == -:
        fivePrimeStartPos = negStrandCigarParser(cigar)
    return fivePrimeStartPos

examples:
softClipper(10M, 0, 10) == 10
softClipper(2S7M, 0, 10) == 8
softClipper(10M, 16, 10) == 19
softClipper(10M5S, 16, 10) == 14

def posStrandCigarParser(cigar(str), position(int)):
''' takes a cigar string and position and calculates the 5' start position, assuming positive strand'''
    figure out if 5' soft clipped
    add soft clip amount to position
    return position

examples:
posStrandCigarParser(10M, 10) == 10
posStrandCigarParser(5S10M, 10) == 5

def negStrandCigarParser(cigar(str), position(int)):
''' takes a cigar string and position and calculates the 5' start position (on the reference genome), assuming negative strand '''
    figure out if 5' soft clipped (end of cigar string)
    calculate genomic "length" of aligned read (handle M, I, D, N, S[r])
    incorporate soft clip amount to position
    return position

negStrandCigarParser(10M, 10) == 19
negStrandCigarParser(10M5S, 10) == 24

def strandParser(strand(int)) -> str:
''' takes flag column 2 in sam file and interprets which strand the read is from '''
    check if read is unmapped. if read is unmapped, strand = unmapped
    if bitwise flag operation says strand is positive, strand = positive
    if bitwise flag operation says strand is negative, strand = negative
    return strand

strandParser(0) == positive
strandParser(4) == unmapped
strandParser(16) == negative
strandParser(256) == secondary

def qnameParser(qname(str)) --> str:
''' takes qname string and parses it to return the umi '''
    return umi

qnameParser(NS500451:154:HWKTMBGXX:1:11101:24260:1121:AAA) == AAA

*** algorithm ***

### first, set empty variables
validUmis = load in umis into a set
uniqueReads = set() # set to hold all unique reads to keep a record of what was seen
dupes = 0 # counter for counting how many duplicate reads
uniques = 0 # counter for counting the unique reads
unmapped = 0 # counter for unmapped reads

### start parsing sam file
open and read in sam file

while true:
### variables necessary for finding a dupe - resets each line of sam file
    cigar = ""
    chrom = 0
    strand = ""
    position = 0
    umi = ""

    line = filehandle.readline()
# always have break condition for a while loop
    if line == "", # when we reach the end of the file...
        BREAK

# start parsing record for info
    if line starts with @:
        continue
    else:
        splitLine = line.strip('/n')
        splitLine = line.split(/t)

# columns of each bit of info that we need
    cigar = column 6
    chrom = rname column 3
    strand = flag col 2
    position = pos col 4
    qname = qname sam col 1

# start translating sam file entries
    whichStrand = strandParser(strand)
    umi = qnameParser(qname)

# decisions about translated entries
    if whichStrand == unmapped, throw out read - don't want unmapped reads
        unmapped += 1
    elif whichStrand != unmapped:
        if umi in validUmis:
            5primeStartPos = softClipper(cigar, whichStrand, position)
            read = set(cigar, chrom, whichStrand, 5primeStartPos, umi)
                if read not in uniqueReads:
                    uniqueReads.append(read)
                    uniques += 1
                    write line to file - build out sam file of unique reads
                else:
                    dupes += 1
    else: continue (whichStrand == unmapped)
