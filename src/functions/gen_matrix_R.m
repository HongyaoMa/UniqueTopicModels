function [R, alpha] = gen_matrix_R(alpha, alpha0)
% gen_MatrixR generates the topic-topic co-occurrence matrix for topic with
% Dirichlet distribution with parameter alpha. alpha 0 is used to normalize
% 

% The Autism NMF Matrix
% Hongyao Ma
% Created:  01/05/2015

if nargin == 1
    alpha0 = sum(alpha);
end

% Renormalize the parameters
alpha = alpha / sum(alpha) * alpha0;

% Find the R matrix: 
% R_{i,i} = alpha_i *(alpha_i + 1) / (alpha_0(alpha_0+1))
% R_{i,j} = alpha_i * alpha_j  (alpha_0 * (alpha_0+1))

R = alpha' * alpha ;
R = R + diag(alpha);
R = R/(alpha0 * (alpha0 + 1));