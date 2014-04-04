function [ A, f_vals ] = solve_lin_alm(M, tau, mu, rho, iterations, tol)
%SOLVE_SUB_GRAD
%   This function solves the following problem
%   
%   min_A    tau * | A |_* 
%   s.t. P.*M = P.*A + P.*E
% 
%   which we convert to the unconstrained problem
% 
%   min_A   tau * | A |_* + < Y, P.*A - P.*M > + mu/2 | P.*A - P.*M |_F^2
% 
%   solved by basic subgradient descent
% 
%   Written by Stephen Tierney

if ~exist('tau', 'var')
    tau = 10;
end

if ~exist('mu', 'var')
    mu = 1;
end

if ~exist('rho', 'var')
    rho = 1;
end

if ~exist('iterations', 'var')
    iterations = 100;
end

if ~exist('tol', 'var')
    tol = 10^-6;
end

f_vals = zeros(iterations, 1);
last_f_value = Inf;

omega = find(M);

P = zeros(size(M));
P(omega) = 1;

A = zeros(size(M));
Y = zeros(size(M));

for k = 1 : iterations
    
    %% Take a step
    partial = mu * (P.*A - (M - 1/mu * Y));
    
    B = A - 1/rho * partial;
    
    [A, s] = nn_prox(B, tau/rho);
    
    %% Step 2, update Y
    Y = Y + mu * (P.*A - M);
    
    %% Check function value
    f_vals(k, 1) = tau * sum(s);
    
    if (abs(last_f_value -  f_vals(k, 1)) <= tol)
        break;
    else
        last_f_value = f_vals(k, 1);
    end
    
end
    
    
end

