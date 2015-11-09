#!/bin/bash

# ---SUMMARY--- Use bedtools intersect to count overlaps between bed12 format read mappings and a set of regions (bed)

# ---TAG--- bed
# ---TAG--- next-gen
# ---TAG--- bioinformatics
# ---TAG--- bedtools

###################################################
#
#       Dave Gerrard, University of Manchester
#       2015
#
###################################################

# usage-message.sh
: ${1?"
        Usage: $0 LIST_FILE REGION_FILE(bed)
        This script expects an input file.
	LIST_FILE:  a text file listing bed12 mapping files to count
	REGION_FILE: a bed file of regions to count into.
	Requires bedtools"}
#  Script exits here if command-line parameter absent,


LIST_FILE=$1
Q1=$2
#Q1=${2-'g18.exons.nonredund.strand.bed'}        # assign optional second command line argument or use default.
#Q2=${3-'UNEX_emb_strand_nonExon.bed'}

Q1_BASE=$( basename $Q1 '.bed' )
#Q2_BASE=$( basename $Q2 '.bed' )

declare -a FILES


#readarray FILES < $LIST_FILE
IFS=$'\n' read -d '' -r -a FILES < $LIST_FILE

for F in ${FILES[@]}; do
        echo $F

       # echo "wget -nc $F"
       # wget -nc -F  $F

        F_BASE=$( basename $F '.bed12' )
        echo "$F_BASE"
	DIR=$( dirname "$F" )
        #echo "gunzip $F_BASE.bed.gz"
        #gunzip $F_BASE.bed.gz

        wc $F >> files.wc.log

        # R approach too slow with memmory problems.
        #echo " Rscript bedOverlapCount.R -q $Q -t $F_BASE.bed -o $F_BASE.bedCount"
        #Rscript bedOverlapCount.R -q $Q -t $F_BASE.bed -o $F_BASE.bedCount

        # using bedtools intersect
        #intersectBed -c -s -a g18.exons.nonredund.strand.bed -b GSM1010938_UCSD.Left_Ventricle.mRNA-Seq.STL001.bed > GSM1010938_UCSD.Left_Ventricle.mRNA-Seq.STL001.intersect.bedCount
        echo "bedtools intersect -c -s  -a $Q1 -b $F  > $DIR/$F_BASE.$Q1_BASE.bedCount"
        bedtools intersect -c -s -split -a $Q1 -b $F  > $DIR/$F_BASE.$Q1_BASE.bedCount

        #echo "bedtools intersect -c -s  -a $Q2 -b $F  > $DIR/$F_BASE.$Q2_BASE.bedCount"
        #bedtools intersect -c -s -split -a $Q2 -b $F  > $DIR/$F_BASE.$Q2_BASE.bedCount

        #echo "rm $F_BASE.bed"
        #rm $F_BASE.bed
done



