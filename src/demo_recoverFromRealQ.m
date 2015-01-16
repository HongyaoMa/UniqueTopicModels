% Generate real Q matrices and recover topic distributions from them

% The Autism NMF Project
% Hongyao Ma
% Created:  01/10/2015
% Modified: 01/12/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

%% Parameters

% A = [1 1 1 0 0 0 0 1 1; ...
%        0 1 1 1 1 1 0 0 0; ...
%        0 0 0 0 1 1 1 1 1]';
%    
% % Sizes
% V = size(A, 1);
% k = size(A, 2);

% % Normalize 
% A = A ./ repmat(sum(A), V, 1);

V = 1000;
k = 5;

% Probability of sets of anchors
p_anchor = [0.001, 0.01, 0.1];
p_anchor = sort(p_anchor);
worst_anchors = 1:k;
best_anchors = (length(p_anchor)-1)*k+1:(length(p_anchor))*k

% Topic distribution matrix A
A = gen_matrix_A(V, k, p_anchor);

% The Dirichlet Parameters
alpha = randi(10, 1, k);
% alpha = [0.1, 0.1, 0.2, 0.3, 0.3];
[R, alpha] = gen_matrix_R(alpha, 0.3);

% The word-word co-occurrence matrix
Q = A * R * A';

%%
% Topics
n_Doc = 2000;
topics = drchrnd(alpha, n_Doc);
R_empirical = topics' * topics / n_Doc;

% Documents
l_Doc = 50;
x = gen_Docs(topics, A, l_Doc);

%% Analyze the matrix Q

% Word-Word Co-occurrance
[Q_emp, Q_bar_emp] = gen_matrix_Q(x, 0);

% Normalized Squared Error
err_Q = norm(Q - Q_emp, 'fro') / norm(Q, 'fro')

% % Rank
% rank_Qemp = rank(Q_emp)
% 
% % Singular Vaules
% [~, S, ~] = svd(Q_emp);
% figure;
% stem(diag(S))

%% Recover the topics

% Find the anchors
candidates = 1:V;
[anchor_inds, anchors] = find_anchors(Q_emp, candidates, k, 0);
% [anchor_inds_conv, anchors] = find_anchors_conv(Q_emp, candidates, k, 0);

% Sort the anchors
anchor_inds = sort(anchor_inds)
anchor_modk = mod(anchor_inds, k);
anchor_modk(anchor_modk == 0) = k
[~, ind_sort] = sort(anchor_modk);
anchor_inds = anchor_inds(ind_sort);

% Recover the matrix A
[A_rec, R_rec] = recoverKL(Q_emp, anchor_inds);
% [A_rec, R_rec] = recover(Q_emp, anchor_inds)

err_A = norm(A - A_rec, 'fro') / norm(A, 'fro')

%% Recover with better/worst set of anchors
[A_rec, R_rec] = recoverKL(Q_emp, best_anchors);
err_A_best = norm(A - A_rec, 'fro') / norm(A, 'fro')

[A_rec, R_rec] = recoverKL(Q_emp, worst_anchors);
err_A_worst = norm(A - A_rec, 'fro') / norm(A, 'fro')


%% With Recover
[A_rec, R_rec] = recover(Q_emp, best_anchors);
err_A_best = norm(A - A_rec, 'fro') / norm(A, 'fro')

[A_rec, R_rec] = recover(Q_emp, worst_anchors);
err_A_worst = norm(A - A_rec, 'fro') / norm(A, 'fro')