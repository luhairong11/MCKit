function [ A, f_vals ] = solve_e_lin_extalm(M, tau, lambda, rho, iterations, tol)
%SOLVE_EXT_GRAD
%   This function solves the following problem
%   
%   min_A    tau * | A |_* + lambda/2 | P.*E |_F^2 
%   s.t. P.*M = P.*A + P.*E
% 
%   which we convert to the unconstrained problem
% 
%   min_A   tau * | A |_* + lambda/2 | P.*A - P.*M |_F^2
% 
%   solved by extended gradient descent
% 
%   Written by Stephen Tierney

if ~exist('tau', 'var')
    tau = 10;
end

if ~exist('lambda', 'var')
    lambda = 1;
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

gamma = 1.1;

f_vals = zeros(iterations, 1);
last_f_value = Inf;

omega = find(M);

P = zeros(size(M));
P(omega) = 1;

A = zeros(size(M));

for k = 1 : iterations
    
    searching = true;
    while( searching )
        
        partial = lambda * (P.*A - M);
        V = A - 1/rho * partial;
        [A, s] = nn_prox(V, tau/rho);
        
        f_vals(k, 1) = tau * sum(s) + lambda/2 * norm(P.*A - M, 'fro')^2;
        
        approx = tau * sum(s) + rho/2 * norm(A - V, 'fro')^2;
        
        if ( f_vals(k, 1) > approx)
            rho = gamma * rho;
        else
            searching = false;
        end
        
    end
    
    if (abs(last_f_value -  f_vals(k, 1)) <= tol)
        break;
    else
        last_f_value = f_vals(k, 1);
    end
    
end

end

