#!/bin/bash

# ---SUMMARY--- Use samtools to divide a bam file into FIRST (Watson) and SECOND (Crick) chromosomal strands.

# ---TAG--- bam
# ---TAG--- bioinformatics
# ---TAG--- next-gen
# ---TAG--- samtools

###################################################
#
#       Dave Gerrard, University of Manchester
#       2015
#
###################################################

# usage-message.sh
: ${1?"
        Usage: $0 FILE 
        Split a .bam file into FIRST and SECOND strand parts.
	Requires samtools.  Output into same dir as input."}
#  Script exits here if command-line parameter absent,


F=$1
#S=${2-'defaultValue'}        # assign optional second command line argument or use default.

BASE=$(basename $F '.bam')
OUT_DIR=$(dirname $F)

FIRST_FILE=$BASE.FIRST.sam
FIRST_SORTED=$BASE.FIRST.sorted.sam
# Allow through headers, take:  flag 80  (64 + 16) = first of pair AND rev-comp
#				flag 128 (2nd of pair)  AND NOT 16 (NOT rev-comp)
samtools view -h $F | gawk 'BEGIN{OFS="\t";}{if(($1 ~ /^@/) || (and($2,80) == 80) || ((and($2,128) == 128) && (and($2,16) != 16))) print $0}' > $OUT_DIR/$FIRST_FILE
#samtools view -H  $F > $OUT_DIR/$FIRST_FILE	# sam file requires header.
#samtools view -f 80 $F >> $OUT_DIR/$FIRST_FILE
#samtools view -f 128 -F 16 $F >> $OUT_DIR/$FIRST_FILE
samtools view -bhS $OUT_DIR/$FIRST_FILE -o $OUT_DIR/$BASE.FIRST.bam
#samtools sort -o $OUT_DIR/$BASE.FIRST.sorted.bam $OUT_DIR/$BASE.FIRST.bam
rm $OUT_DIR/$FIRST_FILE

SECOND_FILE=$BASE.SECOND.sam
# Allow through headers, take:	flag 144 (128 + 16) = second of pair AND rev-comp
#				flag 64 (1st of pair) AND NOT 16 (NOT rev-comp)
samtools view -h $F | gawk '{if(($1 ~ /^@/) || (and($2,144) == 144) || ((and($2,64) == 64) && (and($2,16) != 16))) print $0}' > $OUT_DIR/$SECOND_FILE
#samtools view -H  $F > $OUT_DIR/$SECOND_FILE     # sam file requires header.
#samtools view -f 144 $F >> $OUT_DIR/$SECOND_FILE
#samtools view -f 64 -F 16  $F >> $OUT_DIR/$SECOND_FILE
samtools view -bhS $OUT_DIR/$SECOND_FILE -o $OUT_DIR/$BASE.SECOND.bam
#samtools sort -o $OUT_DIR/$BASE.SECOND.sorted.bam $OUT_DIR/$BASE.SECOND.bam
rm $OUT_DIR/$SECOND_FILE

















## utilites

#for Q in $Q_LIST; do
#done

#BASE=`basename $F '.wig.gz'`
#


#qsub $QSUB_PARAMS -N $BASE perl $SCRIPT_DIR/GetBinsFromBed.pl $DATA_DIR/chr11.1000.bed $BW_DIR/$BASE/$BASE.bw $BASE.chr11.1000.bed.counts
#cd ..
