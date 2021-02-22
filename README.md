# dlrcode
Discrete Lehmann representation codes

This repository contains codes to build and work with the discrete
Lehmann representation for imaginary time Green's functions.



Compile and clean instructions

(1) Simply run ./compile. This will first compile the id_dist library,
which contains subroutines to do pivoted QR and interpolative
decompositions. Make sure all of the tests finish correctly. Note that
on some machines, modifications may need to
be made to the Makefile in id_dist, according to the documentation of
that package. After, this will compile all other code, put object files
into the bin folder, and copy the id_dist static library id_lib.a to the
lib-static folder.

(2) To clean and start from scratch, run ./clean. This will first run
make clean in the id_dist library, and then it will delete all object
files from the bin folder, and the static library id_lib.a from the
lib-static folder.

(3) To run tests, go to the test folder. Each .f90 test program has some
parameters that can be adjusted at the top. It also has a compile script
associated with it called, for example, cmp_testcode, for the test
program testcode.f90. Running the compile script generates an
executable, test, for that test program.



Notes for the user

(1) For the time being, look through the test codes to see how
subroutines are used, along with the documentation of subroutines in the
source files in the src folder for explanations of variables and
subroutines.

(2) Throughout the code, imaginary time (tau) grid points that are in
(0.5,1) are stored in a peculiar manner. Namely, suppose tau in (0.5,1).
Then instead of storing t, we store the number tau' = tau-1. In other words,
we store the negative distance of tau to 1, rather than tau itself. To
recover the correct tau points (for example from the imaginary time DLR
grid dlrit), simply leave any points in [0,1/2] U {1} alone, and add 1 to any
points in (1/2,1).

The reason for this has to do with maintaining full relative accuracy in
floating point arithmetic. To evaluate the kernel K(tau,omega), we
sometimes need to compute the value 1-tau. In particular, we sometimes
need to compute 1-tau for tau very close to 1. If we work with tau
directly, there is an immediate loss of accuracy due to catastrophic
cancellation. If we instead compute tau' to full relative accuracy and
work with it directly rather than with tau, for example by exploiting
symmetries of K(tau,omega) to avoid evaluating 1-tau, we maintain full
relative accuracy. Note, in particular, that tau' cannot be computed by
first computing tau and then subtracting 1; this would defeat the
purpose and lose the digits the use of tau' is intended to maintain.

To avoid confusion, users should therefore use provided subroutines,
which hide this technical complexity, to carry out all operations. In a
situation, for example, in which the user wants to provide a point tau
in (0,1/2) at which to evaluate a DLR, there are two options: (i)
compute tau' to full relative accuracy, and provide this according to
the instructions in the relevant subroutines, thereby maintaining full
relative accuracy in calculations, or (ii) simply convert the provided
tau point to tau' by subtracting 1, or preferably, process all input
points through the subroutine xxxxxx, which will carry out the
conversion for the user. (ii) is the simpler option, but may result in a
mild loss of accuracy in calculations.
