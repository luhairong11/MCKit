function [ H, h_omega ] = get_mean_cols(A, n_cols)
%GET_MEAN_COLS Summary of this function goes here
%   Detailed explanation goes here

empty_cols = randsample(size(A, 2), n_cols);

H_ones = zeros(size(A));
H_ones(:, empty_cols) = 1;
h_omega = find(H_ones);

H = zeros(size(A));
H(h_omega) = A(h_omega);

end

