n1 = 150;
n2 = 300;
r = 10;

randn('state',2009);
rand('state',2009);

A = randn(n1,r) * randn(r,n2);

df = r*(n1+n2-r);

oversampling = 5;
m = min(5*df,round(.99*n1*n2));

omega = randsample(n1*n2,m);

M = zeros(size(A));
M(omega) = A(omega);

[X, f_vals] = solve_e_lin_accalm(M);