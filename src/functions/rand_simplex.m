function x = rand_simplex(d, N)
% rand_simplex generates random points on unit simplex
%
% INPUT
%   d: Dimension of the vectors
%   N: Number of vectors
%
% OUTPUT
%   x: d by N matrix with generated vectors on the unit simplex

% The Autism NMF Project 
% Hongyao Ma
% Created:  10/13/2014
% Modified: 10/14/2014

x = rand(d, N);
colsum = repmat(sum(x), d, 1);
x = x ./ colsum;

end
