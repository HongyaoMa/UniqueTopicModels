% Testing the Topic Recovery Algorithms with correct Q matrices
% Demo or recoverL2
% Demo of Recover

% The Autism NMF Project
% Hongyao Ma
% Created:  12/17/2014
% Modified: 01/01/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

% Load the generated data
% load('syn_Q_1');
% load('syn_Q_tricky4');
% load('syn_Q_50_5_005');
load('syn_Q_500_10_001');
% load('syn_Q_500_20_0001');
% load('syn_Q_500_20_000001_randR');
% load('syn_Q_9_3_001_randR');

%% Find the anchors
candidates = 1:V;
n_anchor = k;
[anchor_inds, anchors] = find_anchors_conv(Q, candidates, n_anchor, 0);
anchor_inds = sort(anchor_inds)
% anchor_inds = [1:9]

%% recoverL2
A
anchor_inds = [1:9]
k = 9;

[A_rec, R] = recoverL2(Q, anchor_inds)

error_KL = norm(A(:, 1:9) - A_rec, 'fro') / norm(A(:, 1:9), 'fro')


%% Recover
[A_rec, R] = recover(Q, anchor_inds)
error_Rec = norm(A(:, 1:9) - A_rec, 'fro') / norm(A(:, 1:9), 'fro')




