function [M0,M1,A,Y,B] = test_matrix(m,n)
%test_matrix Approximate correlation matrix with positive definite leading block.
%   [M0,M1,A,Y,B] = test_matrix(m,n) constructs an m-by-m random correlation
%   matrix A, an n-by-n symmetric random matrix B with unit diagonal, and
%   an m-by-n random Y, such that M0 = [A Y; Y' B] is symmetric indefinite
%   and M1 = blkdiag(A,eye(m)).

A = gallery('randcorr',m);

M1 = blkdiag(A,eye(n));

while 1

    B = -1 + 2.*rand(n,n);    % Values in [-1,1].
    B = triu(B) + triu(B,1)'; % Symmetrize.
    B(1:n+1:n^2) = 1;         % Unit diagonal.

    Y = -1 + 2.*rand(m,n);

    M0 = [A Y; Y' B];

    if min(eig(M0)) < 0, return, end

end
