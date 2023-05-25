%% Least squares example
%
m = 16; n=8;
A = randn(m,n);
b = randn(m,1);

bnds = randn(n,2);
l = min( bnds, [], 2 );
u = max( bnds, [], 2 );

cvx_begin
    variable x(n)   % creates variable x of dimension n
    minimize ( norm(A*x-b) )
    subject to
        l <= x <= u
cvx_end


