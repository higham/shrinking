`shrinking` - MATLAB Codes for Restoring Definiteness via Shrinking
==========

About
-----

`shrinking` is a collection of MATLAB functions for repairing invalid
(indefinite) covariance and correlation matrices, based on the paper

N. J. Higham, N. Strabić, and V. Šego, "[Restoring Definiteness via Shrinking, with an
Application to Correlation Matrices with a Fixed
Block](http://eprints.ma.man.ac.uk/2191/)", MIMS EPrint
2014.54, Manchester Institute for Mathematical
  Sciences, The University of Manchester, UK, November 2014.

The main functions are

* `shrink_bisection`: a bisection algorithm,

* `shrink_bisection_fb`: a bisection algorithm optimized for a block
  2-by-2 matrix with a positive definite (1,1) block A and a 
  block-diagonal target matrix with diagonal blocks A and the identity.

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

* `test_matrix.m`: used by the previous two M-files to generate a test matrix.


Requirements
-------------

The functions `shrink_newton`, `shrink_gep` and `shrink_gep_fb` require
the NAG Toolbox for MATLAB.

The codes have been developed under MATLAB 2014a and 2014b.

License
-------

See `license.txt` for licensing information.

