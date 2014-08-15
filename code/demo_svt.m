paths = [genpath('common'), genpath('standard')];
addpath(paths);

rng(1);

n = 150;
n_blocks = 3;

A = get_block_diag(n, n_blocks);

df = n_blocks*(n + n - n_blocks);

oversampling = 5;
m = min(5*df,round(.99*n*n));

omega = randsample(n*n, m);

M = zeros(size(A));
M(omega) = A(omega);

[X, f_val, stop_val] = solve_svt(M, omega);

rmpath(paths);