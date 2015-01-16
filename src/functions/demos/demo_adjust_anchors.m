% Demo of adjust_anchors

% The Autism NMF Project
% Hongyao Ma
% 01/15/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

A = [eye(3)*1, eye(3)*3, eye(3)*9]';

% Sizes
V = size(A, 1);
k = size(A, 2);

% Normalize
A = A ./ repmat(sum(A), V, 1);

% Topics
alpha = ones(1, k);
R = gen_matrix_R(alpha, 1);

% 2nd Moment
Q = A * R * A';

% Find the anchors
candidates = 1:V;
[anchor_inds, ~] = find_anchors(Q, candidates, k, 0);
anchor_inds

% Adjust the anchors
anchor_inds_adjusted = adjust_anchors(Q, anchor_inds, 0.01)
