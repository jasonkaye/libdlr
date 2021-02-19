# dlrcode
Discrete Lehmann representation codes

This repository contains codes to build and work with the discrete
Lehmann representation for imaginary time Green's functions.



Compile instructions

(1) You first need to compile the id_dist library, which contains
subroutines to do pivoted QR and build interpolative decompositions.
Enter the id_dist directory and doing make.  Make sure all of the tests
finish correctly. Note that on some machines, modifications may need to
be made to the Makefile in id_dist, according to the documentation of
that package.

(2) Do ./compile. This will compile all necessary code, put object files
into the bin folder, and copy the id_dist static library id_lib.a to the
lib-static folder. Ignore all pesky warnings.
