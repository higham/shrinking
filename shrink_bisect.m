function alpha = shrink_bisect(M0,M1,tol)
%shrink_bisect   Shrinking by bisection method.
%   alpha = shrink_bisect(M0,M1,tol) computes the smallest t in [0,1] such
%   that S(t) = t*M1 + (1-t)*M0 is positive semidefinite, where M0 is
%   symmetric indefinite and M1 is symmetric positive definite.
%   tol is a convergence tolerance; default: tol = 1e-4.
%   The purpose of shrinking is to replace the indefinite symmetric 
%   matrix M0 by the positive semidefinite matrix S(alpha).  
%   This function uses the bisection method with Cholesky factorization.
%
%   Note: M0 and M1 are not checked for validity.

if nargin < 3, tol = 1e-4; end

xl = 0;
xr = 1; 
error = xr - xl;
n_iter = 0;

while error > tol
    
   n_iter = n_iter + 1;

   xm = 0.5*(xl + xr);
   S = xm*M1 + (1-xm)*M0;

   [~,p] = chol(S);

   if p  
      xl = xm;  % S is indefinite.
   else
      xr = xm;
   end

   error = xr - xl;

end

alpha = xr;
