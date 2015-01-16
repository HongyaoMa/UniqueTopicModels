function [A, R] = recover(Q, anchor_inds)
% RECOVER finds the topic distribution matrix A and the topic-topic
% co-occurrence matrix R from the word-word co-occurrence matrix Q and a
% given set of anchors
% The algorithm used is the original Recover algorithm 
% Reference: Arora, S., Ge, R., & Moitra, A. (2012). Learning Topic Models 
% -- Going beyond SVD. 2012 IEEE 53rd Annual Symposium on Foundations of 
% Computer Science, 1–10. doi:10.1109/FOCS.2012.49

% The Autism-NMF Project
% Hongyao Ma
% Created:  01/01/2015
% Modified: 01/01/2015

%% The matrix Q_bar
% Q_rowsum = sum(Q, 2);
% Q_bar = Q ./ repmat(Q_rowsum, 1, size(Q, 2));

% Rows corresponding to anchors
Qs = Q(anchor_inds, :);                 % k by V matrix
Qss = Qs(:, anchor_inds);

%% The Recover
Qs_rowsum = sum(Qs, 2);
z = Qss \ Qs_rowsum;
Dinv = diag(z);

A = ((Qss * Dinv) \ Qs)';
R = Dinv * Qss * Dinv;
