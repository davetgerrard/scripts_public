#!/bin/bash

###################################################
#
#       Dave Gerrard, University of Manchester
#       2015
#
###################################################

# usage-message.sh
: ${1?"
        Usage: $0 BED_FILE CHROM_SIZE_FILE N_FIELDS [OUT_DIR]
        This script expects a bed file and chrom size file"}
#  Script exits here if command-line parameter absent,



BED_FILE=$1


#CHROM_SIZE_FILE='/mnt/lustre/scratch/mqbssdgb/k_hanleyRnaSeq/data/genome_table.hg19'
CHROM_SIZE_FILE=$2

N_FIELDS=$3
#FOO=${VARIABLE:-default}
BED_DIR=`dirname $BED_FILE`
BASE_NAME=`basename $BED_FILE .bed`
OUT_DIR=${4:-$BED_DIR}

echo "BED_DIR: $BED_DIR"
echo "N_FIELDS: $N_FIELDS"
echo "BASE_NAME: $BASE_NAME"
echo "OUT_DIR: $OUT_DIR"

tail -n +2 $BED_FILE > $OUT_DIR/$BASE_NAME.bed.temp1




#bedSort $ANALYSIS_DIR/$TOPHAT_DIRPREFIX$SAMPLE/junctions.bed.temp $ANALYSIS_DIR/$TOPHAT_DIRPREFIX$SAMPLE/junctions.bed.temp


bedToBigBed -tab -type=bed$N_FIELDS $OUT_DIR/$BASE_NAME.bed.temp1 $CHROM_SIZE_FILE $OUT_DIR/$BASE_NAME.bb

rm $OUT_DIR/$BASE_NAME.bed.temp1

#rm $ANALYSIS_DIR/$TOPHAT_DIRPREFIX$SAMPLE/junctions.bed.temp

