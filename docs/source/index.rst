.. role:: underline
    :class: underline

Welcome to AFMM2Dlib's Documentation!
**************************************

About :math:`\texttt{AFMM2Dlib}`:
==================================

AFMM (Algebraic FMM) is an FMM implementation where the low rank approximations of appropriate matrix sub-blocks are formed in an algebraic fashion.

:math:`\texttt{AFMM2Dlib}` is a library consisting of AFMM for kernels in 2D.

Low-rank approximation of the appropriate blocks is obtained using ACA. The domain is subdivided based on a KDTree. The algorithm has been parallelized using OpenMP.

The code is written in C++ and features an easy-to-use interface, where the user provides the following inputs:

- a ``kernel`` object which abstracts data of the matrix through a member function ``getMatrixEntry(int i, int j)`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix.

- locations of nodes in the domain through an Eigen matrix ``loc``

- the vector ``b`` to be multiplied to the matrix


The current release has the following capabilities:

- MatVecs: Obtains :math:`A x` at a cost of :math:`\mathcal{O}\left(N\right)`

Doc Contents
============
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   benchmarks

Other Links
===========

Learn more about :math:`\texttt{AFMM2Dlib}` by visiting the

* Code Repository: http://github.com/sivaramambikasaran/AFMM2D
* Documentation: http://afmm2dlib.rtfd.io
