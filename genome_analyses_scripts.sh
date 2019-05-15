#pipeline to analyse different assemblies

###############################
# GENOME ASSEMBLY USING ABYSS #
###############################

#!/bin/bash

#np for number of threads
#k for kmer size

#using multiple libraries for one sample
abyss-pe np=8 k=31 name=path/to/output/prefix lib='lib1 lib2 lib3 lib4 lib5' lib1='path/to/lib1/file.sra' lib2='path/to/lib2/file.sra' lib3='path/to/lib3/file.sra' lib4='path/to/lib4/file.sra' lib5='path/to/lib5/file.sra' &> path/to/verbose.log

#OR
#using one library
abyss-pe np=8 k=31 name=path/to/output/prefix in='path/to/reads_1.fastq.gz path/to/reads_2.fastq.gz' &> path/to/verbose.log

################################
# GENOME ASSEMBLY USING SPADES #
################################

#NOTE: SPAdes is very processing and memory heavy -> use one of the bigger linux boxes

#!/bin/bash

#-t for number of threads
#-m for amount of maximum memory in Mb

spades.py --pe1-1 path/to/lib1/reads_1.fastq.gz --pe1-2 path/to/lib1/reads_2.fastq.gz --pe2-1 path/to/lib2/reads_1.fastq.gz --pe2-2 path/to/lib2/reads_2.fastq.gz -o path/to/output -t 14 -m 200

#############################
# GENOME ASSEMBLY USING SGA #
#############################

#! /bin/bash -x

#this example has two libraries for one sample

cd path/to/reads

IN1=lib1/reads_1.fastq.gz
IN2=lib1/reads_2.fastq.gz
IN3=lib2/reads_1.fastq.gz
IN4=lib2/reads_2.fastq.gz


# Parameters
SGA_BIN=sga
BWA_BIN=bwa
SAMTOOLS_BIN=samtools
BAM2DE_BIN=sga-bam2de.pl
ASTAT_BIN=sga-astat.py
DISTANCE_EST=DistanceEst

# Overlap parameter used for the final assembly. This is the only argument
# to the script
OL=75

# The number of threads to use
CPU=8

# To save memory, we index $D reads at a time then merge the indices together
D=4000000

# Correction k-mer value
CK=31

# The minimum k-mer coverage for the filter step. Each 27-mer
# in the reads must be seen at least this many times
COV_FILTER=2

# Overlap parameter used for FM-merge. This value must be no greater than the minimum
# overlap value you wish to try for the assembly step.
MOL=55

# Parameter for the small repeat resolution algorithm
R=10

# The number of pairs required to link two contigs into a scaffold
MIN_PAIRS=5

# The minimum length of contigs to include in a scaffold
MIN_LENGTH=200

# Distance estimate tolerance when resolving scaffold sequences
SCAFFOLD_TOLERANCE=1

# Turn off collapsing bubbles around indels
MAX_GAP_DIFF=0

#
# Dependency checks
#

# Check the required programs are installed and executable
prog_list="$SGA_BIN $BWA_BIN $SAMTOOLS_BIN $BAM2DE_BIN $DISTANCE_EST $ASTAT_BIN"
for prog in $prog_list; do
    hash $prog 2>/dev/null || { echo "Error $prog not found. Please place $prog on your PATH or update the *_BIN variables in this script"; exit 1; }
done 

# Check the files are found
file_list="$IN1 $IN2 $IN3 $IN4"
for input in $file_list; do
    if [ ! -f $input ]; then
        echo "Error input file $input not found"; exit 1;
    fi
done

# First, preprocess the data to remove ambiguous basecalls
$SGA_BIN preprocess --pe-mode 1 -o output.fastq $IN1 $IN2 $IN3 $IN4

#
# Error correction
#
# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
$SGA_BIN index -a ropebwt -t $CPU --no-reverse output.fastq

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
$SGA_BIN correct -k $CK --discard --learn -t $CPU -o reads.ec.k$CK.fastq output.fastq

#
# Contig assembly
#

# Index the corrected data.
$SGA_BIN index -a ropebwt -t $CPU reads.ec.k$CK.fastq

# Remove exact-match duplicates and reads with low-frequency k-mers
$SGA_BIN filter -x $COV_FILTER -t $CPU --homopolymer-check --low-complexity-check reads.ec.k$CK.fastq

# Merge simple, unbranched chains of vertices
$SGA_BIN fm-merge -m $MOL -t $CPU -o merged.k$CK.fa reads.ec.k$CK.filter.pass.fa

# Build an index of the merged sequences
$SGA_BIN index -d 1000000 -t $CPU merged.k$CK.fa

# Remove any substrings that were generated from the merge process
$SGA_BIN rmdup -t $CPU merged.k$CK.fa

# Compute the structure of the string graph
$SGA_BIN overlap -m $MOL -t $CPU merged.k$CK.rmdup.fa

# Perform the contig assembly without bubble popping
$SGA_BIN assemble -m $OL -g $MAX_GAP_DIFF -r $R -o assemble.m$OL merged.k$CK.rmdup.asqg.gz

#
# Scaffolding/Paired end resolution
# 
CTGS=assemble.m$OL-contigs.fa
GRAPH=assemble.m$OL-graph.asqg.gz

# Realign reads to the contigs
sga-align --name output.pe $CTGS $IN1 $IN2 $IN3 $IN4

# Make contig-contig distance estimates
sga-bam2de.pl -n $MIN_PAIRS --prefix libPE output.pe.bam #received error -> maybe output.pe

# Make contig copy number estimates
sga-astat.py -m $MIN_LENGTH output.pe.refsort.bam > libPE.astat

$SGA_BIN scaffold -m $MIN_LENGTH --pe libPE.de -a libPE.astat -o scaffolds.n$MIN_PAIRS.scaf $CTGS
$SGA_BIN scaffold2fasta -m $MIN_LENGTH -a $GRAPH -o scaffolds.n$MIN_PAIRS.fa -d $SCAFFOLD_TOLERANCE --use-overlap --write-unplaced scaffolds.n$MIN_PAIRS.scaf

#####################################
# GENOME ASSEMBLY USING soapdenovo2 #
#####################################

#NOTE for soapdenovo2 analyses you will need a config file that tells the assembler what to do

###########################
# CONFIG FILE soapdenovo2 #
###########################

#this is an example config file with two libraries for one sample

#maximal read length
max_rd_len=100
[LIB]

#average insert size
avg_ins=500

#if sequence needs to be reversed
reverse_seq=0

#in which part(s) the reads are used
asm_flags=3

#use only first 100 bps of each read
rd_len_cutoff=100

#in which order the reads are used while scaffolding
rank=2

# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3

#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32

#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=path/to/reads_1.fastq.gz
q2=path/to/reads_2.fastq.gz

[LIB] #second library
avg_ins=200
reverse_seq=0
asm_flags=3
rd_len_cutoff=100
rank=1
pair_num_cutoff=3
map_len=32
q1=path/to/reads_1.fastq.gz
q2=path/to/reads_2.fastq.gz

##########################
# ACTUAL RUN soapdenovo2 #
##########################

#!/bin/bash

#-K is kmer size
#-p is number of threads
#-R resolve repears by reads

cd path/to/read_files
SOAPdenovo-63mer all -s path/to/file.config -K 31 -R -p 14 -o path/to/output/dir 1>path/to/file.ass.log 2>path/to/file.ass.err

#-t is number of threads
#-l is maximum read length

GapCloser -b path/to/file.config -a path/to/soapdenovo-output.scafSeq -o path/to/output -t 14 -l 100 1>path/to/file.ass.log 2>path/to/file.ass.err

######################################
# GENOME RE-ASSEMBLY USING redundans #
######################################

#!/bin/bash

cd /home/redundans/src/redundans

./redundans.py -v -i path/to/reads*.fastq.gz -f path/to/scaffolds -o path/to/output -t 8 --log path/to/log

####################################
# GENOME ANNOTATION USING AUGUSTUS #
####################################

#!/bin/bash

augustus --species=coprinus path/to/filtered/scaffolds_fasta > path/to/filtered/name.gff

######################################
# GENE STATISTICS USING GENOME-TOOLS #
######################################

#!/bin/bash

gt gff3 -sort -tidy path/to/augustus.gff | gt stat -genelengthdistri -o path/to/output-genelengthdist.txt
