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
V = 1000;
k = 5;
n_Doc = 2000;
l_Doc = 50;
alpha0 = 0.03;

% Probability of sets of anchors
p_anchor = [0.001, 0.01, 0.1];
p_anchor = sort(p_anchor);
worst_anchors = 1:k;
best_anchors = (length(p_anchor)-1) * k + 1:(length(p_anchor))*k;

% Topic distribution matrix A
A = gen_matrix_A(V, k, p_anchor);

% The Dirichlet Parameters
alpha = gen_alpha(alpha0, k, 'random');
R = gen_matrix_R(alpha);

% The word-word co-occurrence matrix
Q = A * R * A';

%%
% Topics
topics = drchrnd(alpha, n_Doc);
R_empirical = topics' * topics / n_Doc;

% Documents
x = gen_Docs(topics, A, l_Doc);

%% Analyze the matrix Q

% Word-Word Co-occurrance
[Q_emp, Q_bar_emp] = gen_matrix_Q(x, 0);

% Normalized Squared Error
err_Q = norm(Q - Q_emp, 'fro') / norm(Q, 'fro')


%% Recover the topics

% Find the anchors
candidates = 1:V;
[anchor_inds, anchors] = find_anchors(Q_emp, candidates, k, 0);
% [anchor_inds_conv, anchors] = find_anchors_conv(Q_emp, candidates, k, 0);

% Reorder the anchors
anchor_inds_n = reorder_anchors(anchor_inds)
        
% Recover the matrix A
[A_rec, R_rec] = recoverL2(Q_emp, anchor_inds);
% [A_rec, R_rec] = recover(Q_emp, anchor_inds)

err_A = norm(A - A_rec, 'fro') / norm(A, 'fro')

%% Recover with better/worst set of anchors
[A_rec, R_rec] = recoverL2(Q_emp, best_anchors);
err_A_best = norm(A - A_rec, 'fro') / norm(A, 'fro')

[A_rec, R_rec] = recoverL2(Q_emp, worst_anchors);
err_A_worst = norm(A - A_rec, 'fro') / norm(A, 'fro')


%% With Recover
[A_rec, R_rec] = recover(Q_emp, best_anchors);
err_A_best = norm(A - A_rec, 'fro') / norm(A, 'fro')

[A_rec, R_rec] = recover(Q_emp, worst_anchors);
err_A_worst = norm(A - A_rec, 'fro') / norm(A, 'fro')