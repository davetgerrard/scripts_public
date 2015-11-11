#!/bin/bash

# ---SUMMARY--- Use samtools and awk to remove sequence information from a bam file.

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
        Usage: $0 FILE [OUT_FILE]
        Replace SEQ [field 10] and QUAL [field 11] of bam file with '*'.
	Retain headers
	Optional: OUT_FILE for final file name or will add suffix ANON.bam
	Requires samtools.  Output into same dir as input."}
#  Script exits here if command-line parameter absent,


F=$1


BASE=$(basename $F '.bam')

IN_DIR=$(dirname $F)
ALT_OUT=$IN_DIR/$BASE.ANON.bam
OUT_FILE=${2-$ALT_OUT}
OUT_DIR=$(dirname $OUT_FILE)


TEMP_FILE=$OUT_DIR/$BASE.TMP.sam
echo $ALT_OUT
echo $OUT_FILE
echo $TEMP_FILE

samtools view -h $F | awk 'BEGIN{OFS="\t";}{if($1 ~ /^@/) {print $0} else{$10="*" ; $11="*" ; print $0}}' > $TEMP_FILE
samtools view -bhS $TEMP_FILE -o $OUT_FILE
rm $TEMP_FILE
















## utilites

#for Q in $Q_LIST; do
#done

#BASE=`basename $F '.wig.gz'`
#


#qsub $QSUB_PARAMS -N $BASE perl $SCRIPT_DIR/GetBinsFromBed.pl $DATA_DIR/chr11.1000.bed $BW_DIR/$BASE/$BASE.bw $BASE.chr11.1000.bed.counts
#cd ..
