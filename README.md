

This is a fork of the milc_qcd repository. See https://github.com/milc-qcd/milc_qcd

NOTE: Make sure you're on the *develop* branch:

	git checkout develop

The purpose of this fork is simply to develop some small applications that I will use to compare MILC vs. PETSc-Grid.

I have several basic applications which do a simple multiplication of the Dirac operator with a fermion vector. In some cases, I started by copying the ks_spectrum application and making minimal modifications. As such, when compiling those test applications, a lot of unnecessary additional stuff is compiled, but that is fine.

In each case, the output is the vector obtained by mutiplying the Dirac operator by a vector of ones. The different Dirac operators are as follows:

1. diracmatvec: The HISQ Dirac operator
2. diracmatvec_noksphases: The HISQ Dirac operator but without the staggered phases
3. diracmatvec_nofatlinks_noksphases: The HISQ Dirac operator with the Naik term, but with no fat links (i.e. the X and W links are replaced by U links), and no staggered phases
4. diracmatvec_naive_noksphases: The HISQ Dirac operator but without the Naik term, the fat X links, or the staggered phases


To compile, enter the application directory and then

	make diracmatvec
	
To run the test with the 4^4 sample lattice, do:

	./diracmatvec diracmatvec_4444.in

To run the test with the 16^3x32 sample lattice, do:

	./diracmatvec diracmatvec_16161632.in


Note that MILC stores the field by parity, so when printing out a field, the result will be all sites of one parity followed by all sites of the other parity. Parity refers to checkerboard parity i.e. each site is black or white (even or odd). I have added some scripts in diracmatvec_util/ to do basic cleaning, conversion, and comparison between MILC and PETSc solution vectors.




















