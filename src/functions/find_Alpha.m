function [alpha, alpha0] = find_Alpha(R)
% Find the Dirichlet Parameters alpha from topic co-occurrence matrix R

% The Autism NMF Project
% Hongyao Ma
% Created:  01/05/2015
% Modified: 01/05/2015

alpha = sum(R,2);
alpha0 = (1 - R(1,1)/ alpha(1))/ (R(1,1)/ alpha(1) - alpha(1));
alpha = alpha / sum(alpha) * alpha0;



