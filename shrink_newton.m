function alpha = shrink_newton(M0,M1,tolN,tolB)
%shrink_newton   Shrinking by Newton's method.
%   alpha = shrink_newton(M0,M1,tolN,tolB) uses Newton's method to compute
%   the smallest t in [0,1] such that S(t) = t*M1 + (1-t)*M0 is positive
%   semidefinite, where M0 is symmetric indefinite and M1 is symmetric
%   positive definite.
%   tolN is the convergence tolerance for Newton's method; default: tol = 1e-4.
%   The purpose of shrinking is to replace the indefinite symmetric 
%   matrix M0 by the positive semidefinite matrix S(alpha).  
%   In each step, the eigenvector corresponding to the smallest eigenvalue
%   of a symmetric matrix is computed by tridiagonalization followed by
%   bisection.
%   tolB is the convergence tolerance for bisection, with
%   default tolB = 0: see the documentation for nag_lapack_dstebz.
%
%   Note:
%   - M0 and M1 are not checked for validity.
%   - The 64-bit NAG toolbox for MATLAB is required.

if nargin < 4, tolB = 0; end
if nargin < 3, tolN = 1e-4; end 
    
maxit = 500;
x0 = 0; % Starting point

M = M0 - M1;

for n = 0:maxit
    S = x0*M1 + (1-x0)*M0;
    x_min = eigv4newt(S,tolB); % Eigenvector for the smallest eigenvalue
    x1 = x_min'*M0*x_min/(x_min'*M*x_min);
    if abs(x0-x1) <= tolN, alpha = x1; return, end
    x0 = x1;
    if n == maxit 
          error('Not converged in %2.0f iterations', maxit)
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = eigv4newt(A,tolB)
%eigv4newt      Compute the extremal eigenvector for the Newton method.
%   x = eigv4newt(A) computes the eigenvector corresponding to the smallest
%   eigenvalue of the matrix A by tridiagonalization followed by bisection.
%   tolB is the convergence tolerance for bisection; default: tol = 0, as
%   in NAG.

% Find the smallest eigenvalue of A.

% Find tridiagonal T = Q^TAQ, stored in vectors d (diagonal)
% and e (super/subdiagonal).
% nag_lapack_dsytrd = f08fe.
uplo = 'L'; % Which part of A is used to store T: lower triangle here.
[AOut, d, e, tau, info] = nag_lapack_dsytrd(uplo, A);

if info, error('nag_lapack_dsytrd failed.\n'), end

% Find the smallest eigenvalue of T.
range = 'I';
order = 'B';
vl = 0;
vu = 0;
il = int64(1);
iu = il;
abstol = tolB;
% il=iu means we want the smallest eigenvalue
% order='B' so we can use f08jk for eigenvectors
% nag_lapack_dstebz = f08jj
[m, ~, w, iblock, isplit, info] = ...
    nag_lapack_dstebz(range, order, vl, vu, il, iu, abstol, d, e);
% The eigenvalue is in w(1)

if info, error('nag_lapack_dstebz failed.\n'), end

% Find the corresponding eigenvector of T
% nag_lapack_dstein = f08jk
[z,~, info] = nag_lapack_dstein(d, e, m, w, iblock, isplit);
% Eigenvector is z

if info, error(['nag_lapack_dstein failed\n']), end

% Required eigenvector of A is Qz, Q from f08fe, z from f08jk
% nag_lapack_dormtr = f08fg
[zOut, info] = nag_lapack_dormtr('L', uplo, 'N', AOut, tau, z);
% 'L' since Q is on the left of z, 'N' since we have Q not Q^T

if info, error(['nag_lapack_dormtr failed\n']), end

x = zOut;

end
