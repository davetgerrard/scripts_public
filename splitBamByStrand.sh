#!/bin/bash

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
samtools view -H  $F > $OUT_DIR/$FIRST_FILE	# sam file requires header.
samtools view -f 80 $F >> $OUT_DIR/$FIRST_FILE
samtools view -f 128 -F 16 $F >> $OUT_DIR/$FIRST_FILE
#samtools view $F | grep 'XS:A:+' | grep -v -P '^\t'  >> $OUT_DIR/$FIRST_FILE
#samtools view $F | grep 'XS:A:+'  >> $OUT_DIR/$FIRST_FILE
samtools view -bhS $OUT_DIR/$FIRST_FILE -o $OUT_DIR/$BASE.FIRST.bam
samtools sort -o $OUT_DIR/$BASE.FIRST.sorted.bam $OUT_DIR/$BASE.FIRST.bam
rm $OUT_DIR/$FIRST_FILE

SECOND_FILE=$BASE.SECOND.sam
samtools view -H  $F > $OUT_DIR/$SECOND_FILE     # sam file requires header.
samtools view -f 144 $F >> $OUT_DIR/$SECOND_FILE
samtools view -f 64 -F 16  $F >> $OUT_DIR/$SECOND_FILE
#samtools view $F | grep 'XS:A:-' | grep -v -P '^\t'   >> $OUT_DIR/$SECOND_FILE
#samtools view $F | grep 'XS:A:-'  >> $OUT_DIR/$SECOND_FILE
samtools view -bhS $OUT_DIR/$SECOND_FILE -o $OUT_DIR/$BASE.SECOND.bam
samtools sort -o $OUT_DIR/$BASE.SECOND.sorted.bam $OUT_DIR/$BASE.SECOND.bam
rm $OUT_DIR/$SECOND_FILE

















## utilites

#for Q in $Q_LIST; do
#done

#BASE=`basename $F '.wig.gz'`
#


#qsub $QSUB_PARAMS -N $BASE perl $SCRIPT_DIR/GetBinsFromBed.pl $DATA_DIR/chr11.1000.bed $BW_DIR/$BASE/$BASE.bw $BASE.chr11.1000.bed.counts
#cd ..
