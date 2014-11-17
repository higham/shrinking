function test_shrink(m,n)
%test_shrink   Test function for shrinking codes.
%   test_shrink(m,n) does the following things:
%   - Generates a random symmetric indefinite test matrix
%     M0 = [A Y; Y' B] with unit diagonal and A an m-by-m correlation
%     matrix and a target matrix M1 = blkdiag(A,eye(m))
%   - Plots the smallest eigenvalue function to see where approximately
%     the optimal parameter is.
%   - Runs the shrinking methods to compute it.
%   Default: m = 3, n= 4.

if nargin < 2, n = 4; end
if nargin < 1, m = 3; end

[M0,M1,A,Y,B] = test_matrix(m,n);

% Plot the smallest eigenvalue of S(alpha) = alpha*M1 + (1 - alpha)*M0
% The zero of the function is the optimal shrinking parameter alpha_opt

k = 100;
lambda_min = zeros(1,k+1);

for alpha=0:1/k:1
    S = alpha*M1 + (1 - alpha)*M0;
    t = round(alpha*k) + 1;
    lambda_min(1,t) = min(eig(S));
end

x_alpha = linspace(0,1,k+1);

figure
FS = 14;

plot(x_alpha,lambda_min,'k-')
grid
xL = get(gca,'XLim');
line(xL,[0 0],'Color','k','LineStyle','--');
% set(gca,'YTick',0)
xlabel('\alpha','fontsize', FS,'verticalalignment','middle')
options = {'Interpreter','latex'};
ylabel('$\lambda_{\min}$', options{:}, 'rotation',0,'fontsize', FS, ...
       'verticalalignment','top','horizontalalignment','right')
title('Optimal \alpha is value for which \lambda_{min} = 0')

fprintf('The optimal shrinking parameter:\n')
tol = 1e-6;
alpha_bisect = shrink_bisect(M0,M1,tol)
alpha_bisect_fb = shrink_bisect_fb(A,Y,B,tol)
alpha_newton = shrink_newton(M0,M1,tol,tol)
alpha_gep = shrink_gep(M0,M1,tol)
alpha_gep_fb = shrink_gep_fb(A,Y,B,tol)
