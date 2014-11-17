%SHRINK_WORK

nvals = [5:5:50]

for i = 1:length(nvals)
n = nvals(i)
A = full(gallery('tridiag',n,1,1,1));


% Shrinking

theta = 1 + 1/(2*cos(n*pi/(n+1)))

A_shrink = full(gallery('tridiag',n,1-theta,1,1-theta));
d_shrink = sqrt( (2*(n-1)) ) * theta
d_shrink_actual(i) = norm(A - A_shrink,'fro');
fprintf('d_shrink_actual: %9.2e\n', d_shrink_actual(i))

% NCM

mask = (A ~= 0);
[Xf,iter_fixed(i)] = nearcorr_fixed(A,1e-4,[],inf,[],[],mask,0);
d_NCM_tridiag(i) = norm(A - Xf,'fro');
fprintf('d_NCM_tridiag: %9.2e\n',d_NCM_tridiag(i))

[X,iter(i)] = nearcorr(A,1e-6,[],1000,[],[],0);
d_NCM(i) = norm(A - X,'fro');
fprintf('d_NCM: %9.2e\n', d_NCM(i))

% NCM lower bound.
k = 1:n;
e = 1 + 2*cos(k*pi/(n+1));
ineg = find(e<0);
ineg2 = ceil( 2*(n+1)/3 ):n;  % Formula
% if ~isequal(ineg,ineg2), error('Mismatch!'), end 
NCM_lb(i) = norm(e(ineg));
fprintf('NCM_lb: %9.2e\n', NCM_lb(i)) 

end

for i = 1:length(nvals)
    fprintf('%2.0f & %9.2e & %9.2e (%4.0f) &  %9.2e\\\\\n', ...
             nvals(i), d_shrink_actual(i),  ...
             d_NCM_tridiag(i), iter_fixed(i), d_NCM(i) )
end    