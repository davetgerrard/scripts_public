#!/bin/bash

# ---SUMMARY--- Using an ftp url to download an SRX experiment, download the .sra files, unpack to fastq and concatenate them.

# ---TAG--- sra
# ---TAG--- bioinformatics
# ---TAG--- next-gen
# ---TAG--- fastq
# ---TAG--- ascp
# ---TAG--- sra-toolkit

###################################################
#
#       Dave Gerrard, University of Manchester
#       2015
#
###################################################

# usage-message.sh
: ${1?"
        Usage: $0 SRX_URL 
	Downloads all .sra files, extracts fastq files, concatentats the fastq.
	Requires sra-toolkit, aspera. 
	"}
#  Script exits here if command-line parameter absent,


F=$1
#S=${2-'defaultValue'}        # assign optional second command line argument or use default.

echo "input is $F"
BASE=$(basename $F )
OUT_DIR=$(pwd)

ascp -T -l300M -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp-private.ncbi.nlm.nih.gov:$1 $OUT_DIR
echo "Finished download"

FILES=$(find $OUT_DIR/$BASE/ -type f -name '*.sra' )    #> $OUT_DIR/$BASE.fileList


for S in $FILES; do
	echo "Running fastq-dump on $S"
	fastq-dump --split-files  $S
	rm $S
done

echo "Finished fastq-dump"

mkdir $OUT_DIR/fastq/
for S in $FILES; do
	SRA_BASE=$(basename $S '.sra')
	cat ${SRA_BASE}_1.fastq >> $OUT_DIR/fastq/${BASE}_1.fastq
	rm ${SRA_BASE}_1.fastq
	cat ${SRA_BASE}_2.fastq >> $OUT_DIR/fastq/${BASE}_2.fastq
	rm ${SRA_BASE}_2.fastq
done

echo "Finished concatenating fastq files"

gzip -f $OUT_DIR/fastq/${BASE}_1.fastq
gzip -f $OUT_DIR/fastq/${BASE}_2.fastq

echo "Finished gzip'ing fastq files"

rm -R $BASE/
rmdir $BASE

echo "Removed temp SRA folder"

echo "Finished processing experiment $BASE"

#FILES=$(cat )





