#!/bin/bash

# Compare a Petsc and MILC solution vector for nprocs=1
# Note, may need to edit compare_petsc_vecs.py to set lattice size

# Usage: bash compare_vecs.sh petsc_input_vec milc_input_vec

# First clean the two files
bash clean_petsc_vec.sh $1 petsc_vec.tmp
bash clean_milc_vec.sh $2 milc_vec.tmp

# Then compare them
python compare_petsc_milc_vecs.py petsc_vec.tmp milc_vec.tmp

# Clean up
rm *.tmp














