#!/bin/bash

# Take an output vector from Petsc and clean it up for
# input into Python comparison script


file=$1
outfile=$2

# Extract all lines between the lines containing "Vec Object:" and "Mat Object:"
awk '/Vec Object:/{flag=1; next} /Mat Object:/{flag=0} flag' $file > $outfile

# Remove first line
sed -i '1d' $outfile

# Remove all lines starting with "Process" (only relevant for files produced with nprocs>1)
sed -i '/^Process/d' $outfile

# Remove all spaces
sed -i -e 's/ //g' $outfile

# Replace complex i with j
sed -i -e 's/i/j/g' $outfile

















