#!/bin/bash

# ---SUMMARY--- Use bedtools bamToBed to make a bed12 file but force the strand attribute to a certain value.

# ---TAG--- bam
# ---TAG--- bed
# ---TAG--- bioinformatics
# ---TAG--- bedtools
# ---TAG--- next-gen

###################################################
#
#       Dave Gerrard, University of Manchester
#       2015
#
###################################################

# usage-message.sh
: ${1?"
        Usage: $0 FILE 
        This script expects an input file."}
#  Script exits here if command-line parameter absent,


F=$1
S=${2-'+'}        # assign optional second command line argument or use default.
MIN_QUAL=50

echo "F: $F"

F_BASE=$( basename "$F" .bam )
DIR=$( dirname "$F" )

echo "F_BASE: $F_BASE"
echo "DIR: $DIR"

#bamToBed -i  $F > $DIR/$F_BASE.bed
#gzip $DIR/$F_BASE.bed

# currently prints out only (tophat) uniquely mapped reads.
bamToBed -bed12 -i  $F |  awk -v strand=$S -v minQual=$MIN_QUAL '{if($5 >= minQual) print  $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"strand"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > $DIR/$F_BASE.unique.bed12
#cat $DIR/$F_BASE.bed12 | awk '{if($5 == 50) print $0}' > $DIR/$F_BASE.unique.bed12

#gzip $DIR/$F_BASE.bed


#bamToBed -split -i $F > $DIR/$F_BASE.split.bed
#gzip $DIR/$F_BASE.bed








