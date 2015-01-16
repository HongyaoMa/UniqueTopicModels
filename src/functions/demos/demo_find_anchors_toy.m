% Toy Demo of the find anchors function
% Run find_anchors algorithms just on the basis and plot the results
% Compare the original find_anchor algorithm and find_anchor_conv which
% computes the convex expansions

% The NMF Autism Project
% 12/15/2014

clear;
clc;
close all;

% Add all subfolders of the source and data directory
addpath(genpath('../src'));
addpath(genpath('../data'));

% Loading Data
load('find_anchor_data_5');

% Choose the basis to work with
basis = B5';

%% find_anchors
n_anchor = 3;
candidates = 1:size(basis, 1);
[anchor_indices, anchors] = find_anchors(basis, candidates, n_anchor, 0);


figure;
hold on;
title('FindAnchors');

scatter3(B(1,:), B(2,:), B(3,:), 'ks');
scatter3(basis(:,1), basis(:,2), basis(:,3), 'k*');
scatter3(anchors(:, 1), anchors(:, 2), anchors(:, 3), 'ro');

% scatter3(anchors(2, 1), anchors(2, 2), anchors(2, 3), 'gs');

sum(anchors)

%% find_anchors_conv
n_anchor = 6;
[anchor_indices, anchors] = find_anchors_conv(basis, candidates, n_anchor, 0);

figure;
hold on;
title('FindAnchors-Conv');
scatter3(B(1,:), B(2,:), B(3,:), 'ks');
scatter3(basis(:,1), basis(:,2), basis(:,3), 'k*');

scatter3(anchors(:, 1), anchors(:, 2), anchors(:, 3), 'ro');