function [ X, s ] = nn_prox( Y, tau )
%NN_PROX
%   min lamYda*|X|_* + 1/2*|X - Y|^2

[U, S, V] = svd(Y);

s = diag(S);

ind = find(s <= tau);
s(ind) = 0;

ind = find(s > tau);
s(ind) = s(ind) - tau;

S = diag(s);

if (size(Y,1) ~= size(Y,2))
    rows = size(Y,1);
    cols = size(Y,2);
    S(:,rows+1:cols) = zeros(rows,cols - rows);
end

X = U*S*(V');

end

