function alpha = shrink_gep_fb(A,Y,B,tol)
%shrink_gep_fb   Shrinking by generalized e'value problem; fixed block variant.
%   alpha = shrink_gep_fb(A,Y,B,tol) computes the smallest t in [0,1]
%   such that S(t) = t*M1 + (1-t)*M0 is positive semidefinite, where
%   M0 = [A Y; Y' B] is symmetric indefinite with unit diagonal, A is positive
%   definite and M1 = diag(A,I).
%   The purpose of shrinking is to replace the indefinite symmetric matrix
%   M0 by the positive semidefinite matrix S(alpha); in this variant the
%   (1,1) block A of M0 is fixed.  
%   This function finds mu, the smallest eigenvalue of the symmetric
%   definite pencil M0 - mu*M1 by exploiting the structure and using
%   tridiagonalization followed by bisection, and sets alpha = mu/(mu-1).
%   tol is the convergence tolerance for bisection, with
%   default tol = 0: see the documentation for nag_lapack_dstebz.
%
%   Note:
%   - M0 and M1 are not checked for validity.
%   - The 64-bit NAG toolbox for MATLAB is required.

if nargin < 4, tol = 0; end

m = length(A);
n = length(B);

% Cholesky needed only if M1 is not the identity matrix.
M1 = blkdiag(A,eye(m));
if norm(M1 - eye(size(M1)),1)    
    R = chol(A);
    X = R'\Y;
    C = [eye(m) X; X' B];
else
    C = [A Y; Y' B];
end

% Find the smallest eigenvalue of C.

% Find tridiagonal T = Q^TCQ, stored in vectors d (diagonal)
% and e (super/subdiagonal).
% nag_lapack_dsytrd = f08fe.
uplo = 'L'; % Which part of C is used to store T: lower triangle here.
[~, d, e, ~, info] = nag_lapack_dsytrd(uplo, C);

if info, error('nag_lapack_dsytrd failed.\n'), end

% Find the smallest eigenvalue of T.
range = 'I';
order = 'B';
vl = 0;
vu = 0;
il = int64(1);
iu = il;
abstol = tol;
% il = iu means we want the smallest eigenvalue
% nag_lapack_dstebz = f08jj
[~, ~, w, ~, ~, info] = ...
    nag_lapack_dstebz(range, order, vl, vu, il, iu, abstol, d, e);
% The eigenvalue is in w(1)

if info, error('nag_lapack_dstebz failed.\n'), end

alpha = w(1)/(w(1)-1);