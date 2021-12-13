This is a 2D FMM Code to perform Matrix-Vector products of the form Ax.

It is a completely algebraic version. Compressions are made via ACA.

It is applicable to all FMM-able and symmetric kernels.

The kernel that defines the matrix A is to be defined in function "getMatrixEntry" of "userkernel" class in "kernel.hpp" file.

The vector x is to be defined in function "chargesFunction" of "userkernel" class in "kernel.hpp" file.

"ACA.hpp" file contains the ACA module.

"FMM2DTree.hpp" file contains the algorithm.

Before running make sure Eigen and openmp library paths are specified in Makefile.

The algorithm is built on a KD Tree. It generates uniform tree - each leaf is at level nLevels and the number of particles in boxes at a given level differ by atmost 1.

It takes these inputs at run time: Number of particles in the domain, minimum number of particles in each leaf, tolerance of ACA in powers of 10

Half side length of the square domain is assumed to be 1.

To run it input in terminal:

make -f Makefile2D.mk clean

make -f Makefile2D.mk

The output looks like:

./testFMM2D 102400 36 8

Number of Levels of tree: 5

Number of particles is: 102400

Time taken to create the tree is: 0.001077

Time taken to assemble is: 5.2894

Time taken to do Mat-Vec product is: 0.651081

Performing error calculation in box: 108

err: 9.19835e-08
