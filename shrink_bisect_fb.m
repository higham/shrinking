function alpha = shrink_bisect_fb(A,Y,B,tol)
%shrink_bisect_fb   Shrinking by bisection method; fixed block variant.
%   alpha = shrink_bisect_fb(A,Y,B,tol) computes the smallest t in
%   [0,1] such that S(t) = t*M1 + (1-t)*M0 is positive semidefinite, where
%   M0 = [A Y; Y' B] is symmetric indefinite with unit diagonal, A is
%   positive definite and M1 = diag(A,I).
%   tol is a convergence tolerance; default: tol = 1e-4.
%   The purpose of shrinking is to replace the indefinite symmetric matrix
%   M0 by the positive semidefinite matrix S(alpha); in this variant the
%   (1,1) block A of M0 is fixed.  
%   This function uses the bisection method with Cholesky factorization,
%   optimized to exploit the structure.
%
%   Note: M0 and M1 are not checked for validity.

if nargin < 3, tol = 1e-4; end

xl = 0;
xr = 1;
error = xr - xl;
n_iter = 0;
 
R11 = chol(A);
X = R11'\Y;
X = X'*X;

while error > tol
    
    n_iter = n_iter + 1;
    
    xm = 0.5*(xl + xr);
    S = xm*eye(length(B)) + (1-xm)*B - (1-xm)^2*X;

    [~,p] = chol(S);

    if p 
       xl = xm; % S is indefinite
    else
       xr = xm;
    end

    error = xr - xl;
end

alpha = xr;
