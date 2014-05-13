function [ A, f_vals, stop_vals ] = solve_ialm(M, omega, tau, mu, iterations, tol )
%SOLVE_IALM
%   This function solves the following problem
%   
%   min_A    tau * | A |_*
%   s.t. M = A + E, P.*E = 0
%   
%   E is not actually noise. Rather it is an auxillary variable
%   to enforce the MC constraint. We expect P.*A = P.*M.
% 
%   The method is from
%   "The Augmented Lagrange Multiplier Method for Exact Recovery of Corrupted Low-Rank Matricesn"
%   by Zhouchen Lin, Minming Chen, Yi Ma
% 
%   Written by Stephen Tierney

if ~exist('omega', 'var')
    error('Aborted: no observation set provided.');
end

if ~exist('mu', 'var')
%     mu = 1/norm(M, 2);
    mu = 0.1;
end

if ~exist('tau', 'var')
    tau = 10;
end

if ~exist('iterations', 'var')
    iterations = 100;
end

if ~exist('tol', 'var')
    tol = 10^-6;
end

f_vals = zeros(iterations, 1);
stop_vals = zeros(iterations, 2);

P = zeros(size(M));
P(omega) = 1;

Y = zeros(size(M));
E = zeros(size(M));
R = ones(size(P)) - P;

for k = 1 : iterations

    %% Step 1, solve for A
    
    [A, s] = nn_prox(M - E + (1/mu)*Y, tau/mu);
    
    %% Step 2, update E
    old_E = E;
    E = R.*(M - A + (1/mu)*Y);
    
    %% Step 3, update Y
    
    Y = Y + mu*(M - A - E);
    
    %% Get function value
    f_vals(k, 1) = tau * sum(s);
    
    %% Check stopping criteria
    stop_vals(k, 1) = norm(M - A - E, 'fro') / norm(M, 'fro');
    stop_vals(k, 2) = min(mu, sqrt(mu)) * norm(E - old_E, 'fro') / norm(M, 'fro');
    
    if ( stop_vals(k, 1) <= tol && stop_vals(k, 2) <= tol)
        break;
    end
    
end
    
    
end

