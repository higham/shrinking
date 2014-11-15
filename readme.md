`shrinking` - MATLAB Codes for Restoring Definiteness via Shrinking
===

About
=====

`shrinking` is a collection of MATLAB functions for repairing invalid
(indefinite) covariance and correlation matrices, based on the paper
Higham, Strabić, Šego, "[Restoring Definiteness via Shrinking, with an
Application to Correlation Matrices with a Fixed
Block](http://eprints.ma.man.ac.uk/2191/)"

The module requires `scipy.linalg` and it incorporates the following methods:

* `shrink_bisection`: a bisection algorithm,

* `shrink_bisection_fb`: a bisection algorithm optimized for a block
  2-by-2 matrix with a positive definite (1,1) block A and a 
  block-diagonal target matrix with diagonal blocks A dna the identity.

* `shrink_newton`: a Newton algorithm.

* `shrink_gep`: a generalized eigenvalue-based algorithm.

* `shrink_gep_fb`: a generalized eigenvalue-based algorithm optimized for
  the same fixed block case as `shrink_bisection_fb`.

Other M-files:

* `test_shrink`: tests the above functions on a single test problem and
  plots the function of the underlying optimization problem.

* `test_shrink_timings`: runs timing tests on the shrinking codes as well
  as the NAG code `g02aa/nag_nearest_correlation` for computing the nearest
  correlation matrix.

It is possible to incorporate weights into the target matrix that reflect
the confidence with which individual matrix entries are known. See the
paper above for details on how to do this.


Requirements
=============

The functions `shrink_newton`, `shrink_gep` and `shrink_gep_fb` require
the NAG Toolbox for MATLAB.

The codes have been developed under MATLAB 2014a and 2014b.

License
=======

See `license.txt` for licensing information.

