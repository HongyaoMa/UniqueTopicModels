% Demo for different A and similar Q
% Prof. Finale's Email

% The Autism NMF Project
% Hongyao Ma
% 01/08/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

%% Problem Setup

% Topic Distributions
A1 = [ 1 1 1 1  0 0 0 0  0 0 0 0 ;...
       0 0 0 0  1 1 1 1  0 0 0 0 ;...
       0 0 0 0  0 0 0 0  1 1 1 1 ]';

A2 = [ 1 1 1 1  0 0 0 0  0 0 0 0 ;...
       0 0 0 0  1 1 0 0  1 1 0 0 ;...
       0 0 0 0  0 0 1 1  0 0 1 1 ]';
   
% Topic-Topic Co-occurrance
R = [ 1 0       0 ;...
      0 0.55    0.45 ;...
      0 0.45    0.55 ];
  
% Parameters
V = size(A1, 1);
k = size(A1, 2);

% Normalization
A1 = A1 ./ repmat(sum(A1), V, 1);
A2 = A2 ./ repmat(sum(A2), V, 1);

R = R / sum(sum(R));

%% Q1 + Q2 / 2

% Co-occurrence matrices
Q1 = A1 * R * A1'
Q2 = A2 * R * A2'

% || Q1 - Q2 ||_1
error_l1 = sum(sum(abs(Q1 - Q2)))

% The "noisy" Q
Q = (Q1 + Q2)/2;
sumQ = sum(sum(Q))

% Find the anchors
candidates = 1:V;
rankQ = rank(Q)
[anchor_inds_conv, anchors, rankplus, moreanchors, moreinds] = find_anchors_conv(Q, candidates, rankQ, 0);
[anchor_inds, anchors] = find_anchors_conv(Q, candidates, rankQ, 0);

% The true non-negative rank
rankplus_Q = rankplus

% Sort the anchors
anchor_inds_conv
anchor_inds
moreinds
anchor_inds = sort(anchor_inds);
moreinds = sort(moreinds);

%% Recover the topics

% Recover with 3 anchors
[A_rec, R_rec] = recoverL2(Q, anchor_inds)

[A_rec, R_rec] = recover(Q, anchor_inds)

% Recover with 5 anchors
[A_rec, R_rec] = recoverL2(Q, moreinds)

[A_rec, R_rec] = recover(Q, moreinds)

% % Use find_anchors instead
% [anchor_inds, anchors] = find_anchors(Q, candidates, rankQ, 0);
% [A_rec, R_rec] = recoverL2(Q, anchor_inds)
% [A_rec, R_rec] = recover(Q, anchor_inds)

%% Linear w.r.t. R

alpha = [2 3 5];
R2 = gen_MatrixR(alpha)
Q3 = A1 * R2 * A1';

% Different R
Q = (Q1 + Q3)/2

% Find the anchors
rankQ = rank(Q)
[anchor_inds, anchors, rankplus, moreanchors, moreinds] = find_anchors_conv(Q, candidates, rankQ, 0);

% Sort the anchors
anchor_inds = sort(anchor_inds)
moreinds = sort(moreinds)

% Recover the topic distributions
[A_rec, R_rec] = recoverL2(Q, anchor_inds)
[A_rec, R_rec] = recover(Q, anchor_inds)

% Error of the recovered R
error_R = norm((R + R2)/2 - R_rec, 'fro')

%% A = A1 + A2
% Consider the first case when R is from the Dirichlet distribution

A = (A1 + A2)/2;
Q = A * R2 * A'
sum(sum(Q))

% Find the anchors
rankQ = rank(Q)
[anchor_inds, anchors, rankplus, moreanchors, moreinds] = find_anchors_conv(Q, candidates, rankQ, 0);

% Sort the anchors
anchor_inds = sort(anchor_inds)
moreinds = sort(moreinds)

% Recover the topic distributions
[A_rec, R_rec] = recoverL2(Q, anchor_inds)
[A_rec, R_rec] = recover(Q, anchor_inds)

% Recover with 5 anchors
[A_rec, R_rec] = recoverL2(Q, moreinds)
[A_rec, R_rec] = recover(Q, moreinds)

%% A = A1 + A2
% We have a "bad" R

A = (A1 + A2)/2;
Q = A * R * A'
sum(sum(Q))

% Find the anchors
rankQ = rank(Q)
[anchor_inds, anchors, rankplus, moreanchors, moreinds] = find_anchors_conv(Q, candidates, rankQ, 0);

% Sort the anchors
anchor_inds = sort(anchor_inds)
moreinds = sort(moreinds)

% Recover the topic distributions
[A_rec, R_rec] = recoverL2(Q, anchor_inds)
[A_rec, R_rec] = recover(Q, anchor_inds)

% Recover with 5 anchors
[A_rec, R_rec] = recoverL2(Q, moreinds)
[A_rec, R_rec] = recover(Q, moreinds)

%% Find_anchors
[anchor_inds, anchors] = find_anchors(Q, candidates, rankQ, 0);

anchor_inds

[A_rec, R_rec] = recoverL2(Q, anchor_inds)
[A_rec, R_rec] = recover(Q, anchor_inds)

%% Check the convex hull
Q_bar = Q./ repmat(sum(Q,2), 1, V);
[ind_in, res_norm, coeffs] = inConvHull(Q_bar(anchor_inds, :)', Q_bar')

%% Compare the anchors
[anchor_inds, anchors] = find_anchors(Q, candidates, rankQ, 0);
[anchor_inds_conv, anchors, rankplus, moreanchors, moreinds] = find_anchors_conv(Q, candidates, rankQ, 0);

anchor_inds

anchor_inds_conv

correct_inds = Q_bar(anchor_inds, :)' \ Q_bar'

Qbar1 = Q_bar(1,:);
Qbar5 = Q_bar(5,:);
Qbar7 = Q_bar(7,:);
Qbar1 = Q_bar(11,:);

norm51 = norm(Qbar5 - Qbar1)
norm71 = norm(Qbar7 - Qbar1)