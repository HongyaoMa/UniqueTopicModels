function x_out = colStoc(x_in)
% colStoc converts the input matrix x_in into a column stochastic matrix
% x_out by normalizing each column by its summation

% The Autism NMF 
% Hongyao Ma
% 10/13/2014

colsum = sum(x_in);
x_out = x_in ./ repmat(colsum, size(x_in, 1), 1);