function test_shrink_timings(imax,jmax)
%test_shrink_timings  Timing tests for the shrinking codes and NCM.
%   test_shrink_timings(imax,jmax) does the following:
%    - Generates jmax 2x2-block test matrices with diagonal blocks of
%      order 100*i, for each i from 1 to imax.
%    - For each of the imax groups, outputs the average run time for the
%      shrinking methods and for computing the NCM by NAG code g00aa.
%   Default: imax = 5, jmax = 5.

if nargin < 2, jmax = 5; end
if nargin < 1, imax = 5; end

bisect = zeros(1,imax);
bisect_opt = zeros(1,imax);
newt = zeros(1,imax);
ge = zeros(1,imax);
ge_opt = zeros(1,imax);
ncm = zeros(1,imax);       % Requires NAG toolbox: g02aa,

tol = 1e-6;
tolB = 1e-6;
p = 100;

for i = 1:imax
    fprintf('Group %1.0f of %1.0f\n', i, imax)
    for j = 1:jmax
        fprintf('    Matrix %1.0f of %1.0f\n', j, jmax)
        
        [M0,M1,A,Y,B] = test_matrix(p*i,p*i);
        
        tic
        shrink_bisect(M0,M1,tol);
        t = toc;        
        bisect(i) = bisect(i) + t;
        
        tic
        shrink_bisect_fb(A,Y,B,tol);
        t = toc;        
        bisect_opt(i) = bisect_opt(i) + t;
        
        tic
        shrink_newton(M0,M1,tol,tolB);
        t = toc;
        newt(i) = newt(i) + t;
        
        tic
        shrink_gep_fb(A,Y,B,tolB);
        t = toc;
        ge_opt(i) = ge_opt(i) + t;
        
        tic
        shrink_gep(M0,M1,tolB);
        t = toc;
        ge(i) = ge(i) + t;
        
        tic
        g02aa(M0,'errtol',tol);
        t = toc;
        ncm(i) = ncm(i) + t;        
                
    end
end

bisect = bisect/jmax;
bisect_opt = bisect_opt/jmax;
newt = newt/jmax;
ge = ge/jmax;
ge_opt = ge_opt/jmax;
ncm = ncm/jmax;

fprintf('  (m,n)       bisect    bisect_fb    GEP       GEP_fb     newton       NCM \n')
for i = 1:imax
fprintf('(%d,%d) %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',...
    p*i, p*i, bisect(i), bisect_opt(i), newt(i), ge(i), ge_opt(i), ncm(i))
end