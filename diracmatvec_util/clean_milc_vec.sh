#!/bin/bash

# Take an output vector from MILC and clean it up for
# input into Python comparison script


file=$1
outfile=$2

# Extract all lines between the lines containing "Solution vector:" and "DONE"
awk '/Solution vector:/{flag=1; next} /DONE/{flag=0} flag' $file > $outfile

# Remove all spaces
sed -i -e 's/ //g' $outfile

# Replace complex i with j
sed -i -e 's/i/j/g' $outfile

















