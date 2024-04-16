#! /usr/bin/env python
"""

Leon Hostetler, Nov. 24, 2023

USAGE example: 

	python compare_petsc_milc_vecs.py petsc_vec.out milc_vec.out

"""

def IndexToCoord(index, lat_sizes):

  coord = []
  for d in range(len(lat_sizes)):
    coord.append(index % lat_sizes[d])
    index = int(index / lat_sizes[d])

  return coord


def CoordToIndex(coord, lat_sizes):
  index=0
  mult=1
  for d in range(len(lat_sizes)):
    index += mult*coord[d]
    mult *= lat_sizes[d]

  return index


def GlobalCoorToProcessorCoorLocalCoor(global_coor, global_lat_sizes, procs):

  proc_coor = []
  loc_coor = []

  for d in range(len(global_lat_sizes)):
    tmp = global_lat_sizes[d]/procs[d]
    proc_coor.append(int(global_coor[d]/tmp))
    loc_coor.append(int(global_coor[d]%tmp))
    
  return proc_coor, loc_coor



import numpy as np
import sys

ndim = 4
lattice_sizes = [16,16,16,32]
nums_per_site = 3 # The number of complex numbers associated with each lattice site

nsites = np.prod(lattice_sizes)

# Load vectors
petsc_data = np.genfromtxt(sys.argv[1], comments='#', dtype=complex)
milc_data = np.genfromtxt(sys.argv[2], comments='#', dtype=complex)

# Perform some basic checking
if(len(petsc_data) != nums_per_site*nsites):
  print("ERROR: The length of the vector seems to be incorrect")
if(len(milc_data) != nums_per_site*nsites):
  print("ERROR: The length of the vector seems to be incorrect")
  
milc_vec = milc_data

# Easiest approach will be to convert the Petsc vector to a MILC vector
# By having the first half contain the even parity sites and the last
# half contain the odd parity sites
petsc_even = []
petsc_odd = []
for site in range(nsites):
    cartcoord = IndexToCoord(site, lattice_sizes)
    parity = np.sum(cartcoord)%2 # 0 = even, 1 = odd
    if(parity==0):
        for i in range(nums_per_site):
            petsc_even.append(petsc_data[nums_per_site*site+i])
    else:
        for i in range(nums_per_site):
            petsc_odd.append(petsc_data[nums_per_site*site+i])

petsc_vec = np.concatenate([petsc_even, petsc_odd])

same = True
if (not np.array_equal(petsc_vec,milc_vec)):
    print("WARNING: The vectors are not exactly the same.")
    same=False

diff = petsc_vec-milc_vec
print("The difference vector has the following characteristics:")
print(" L2 norm = ", np.linalg.norm(diff))
print(" Max real part = ", diff.real.max())
print(" Max imaginary part = ", diff.imag.max())
print(" Max absolute value = ", np.abs(diff).max())
