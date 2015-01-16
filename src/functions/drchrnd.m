function r = drchrnd(alpha, n)
% DRCHRND generates n random vectors with distribution Dir(alpha)
% Reference:https://cxwangyi.wordpress.com/2009/03/18/to-generate-random-numbers-from-a-dirichlet-distribution/

% The Autism NMF Project
% Hongyao Ma
% 01/06/2015

% take a sample from a dirichlet distribution
p = length(alpha);
r = gamrnd(repmat(alpha, n, 1), 1, n, p);
r = r ./ repmat(sum(r,2),1,p);