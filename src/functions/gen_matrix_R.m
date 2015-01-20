function R = gen_matrix_R(alpha)
% gen_MatrixR generates the topic-topic co-occurrence matrix for topic with
% Dirichlet distribution with parameter alpha. alpha 0 is used to normalize
% 

% The Autism NMF Matrix
% Hongyao Ma
% Created:  01/05/2015

% Find the R matrix: 
% R_{i,i} = alpha_i *(alpha_i + 1) / (alpha_0(alpha_0+1))
% R_{i,j} = alpha_i * alpha_j  (alpha_0 * (alpha_0+1))

alpha0 = sum(alpha);

R = alpha' * alpha ;
R = R + diag(alpha);
R = R/(alpha0 * (alpha0 + 1));