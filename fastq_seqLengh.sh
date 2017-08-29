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
S=${2-'defaultValue'}        # assign optional second command line argument or use default.


printf "%s\t" $F
zcat $F | awk '{if(NR%4==2) print length($0)}'  | head -n 1




