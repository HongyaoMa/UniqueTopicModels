function x_out = rowStoc(x_in)
% rowStoc converts the input matrix x_in into a row-stochastic matrix by
% normalizing each rwo with its sum

% The Autism NMF Project
% Hongyao Ma
% 10/13/2014

rowsum = sum(x_in, 2);
x_out = x_in ./ repmat(rowsum, 1, size(x_in, 2));