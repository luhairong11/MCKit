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

A_noisey = A + (0.05 * randn(size(A)));

M = zeros(size(A));
M(omega) = A_noisey(omega);

tic; [X, f_vals] = solve_e_lin_ext(M, omega, 3, 1.5); toc;

rmpath(paths);

figure, imagesc(A_noisey);
figure, imagesc(M);
figure, imagesc(X);