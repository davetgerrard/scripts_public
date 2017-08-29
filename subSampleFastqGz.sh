#!/bin/bash

###################################################
#
#       Dave Gerrard, University of Manchester
#       2016
#
###################################################

# usage-message.sh
: ${1?"
        Usage: $0 FILE 
        This script expects an input file."}
#  Script exits here if command-line parameter absent,


F=$1
S=${2-1000000''}        # assign optional second command line argument or use default.




BASE=`basename $F '.fastq.gz'`
DIR=`dirname $F `
zcat  $F  |  awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | shuf | head -n $S | sed 's/\t\t/\n/g' | gzip -9 > $DIR/$BASE.$S.fastq.gz







## utilites

#for Q in $Q_LIST; do
#done

#BASE=`basename $F '.wig.gz'`
#
#QSUB_PARAMS='-V -b y -cwd -q node.q '
# do bin counts in a file specific directory.
#qsub $QSUB_PARAMS -N $BASE perl $SCRIPT_DIR/GetBinsFromBed.pl $DATA_DIR/chr11.1000.bed $BW_DIR/$BASE/$BASE.bw $BASE.chr11.1000.bed.counts
#cd ..
