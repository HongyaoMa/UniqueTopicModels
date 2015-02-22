function [A_rec, R] = recoverL2(Q, anchor_inds, true_anchors)
% recoverL2 finds the topic distribution matrix A and the topic-topic
% co-occurrence matrix R from the word-word co-occurrence matrix Q and a
% given set of anchors
% The algorithm used is the RecoverL2 algorithm 
% Reference: Arora, S., Ge, R., Halpern, Y. Y., Mimno, D., Moitra, A., 
% Sontag, D., … Zhu, M. (2013). A Practical Algorithm for Topic Modeling 
% with Provable Guarantees. In Proceedings of The 30th International 
% Conference on Machine Learning (pp. 280–288). Learning; Data Structures 
% and Algorithms; Machine Learning. Retrieved from 
% http://arxiv.org/abs/1212.4777

% The Autism-NMF Project
% Hongyao Ma
% Created:  01/01/2015
% Modified: 01/01/2015



%% The matrix Q_bar
Q_rowsum = sum(Q, 2);
Q_bar = rowStoc(Q);

% Rows corresponding to anchors
size(Q_bar)
anchor_inds
Qs_bar = Q_bar(anchor_inds, :);                 % k by V matrix

% Recover using true anchors (normalized)
if nargin == 3
    Qs_bar = rowStoc(true_anchors);
end

%% Recover KL
% Find the convex coefficients C
[~, ~, C] = inConvHull(Qs_bar', Q_bar');

C = C';
A_prime = diag(Q_rowsum) * C;

% Normalize the matrix A
pz_rec = sum(A_prime);
A_rec = A_prime ./ repmat(pz_rec, size(Q, 2), 1);

%% Find the matrix R
A_pinv = pinv(A_rec);
R = A_pinv * Q * A_pinv';


%% TODO
% 1. Sometimes the A_rec would have something like nan of inf in it thus
% the pinv at line 38 would generate an error
