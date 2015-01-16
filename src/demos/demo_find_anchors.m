% Demo find_anchors and find_anchors_conv in 3 dimensions
% Plot the residue at each step
% Visualize the anchors found at each step
% find_anchors is unable to deal with the cases where the rank is smaller
% than the non-negative rakn however find_anchors_conv can keep looking
% for more anchors than the rank of the data

% The Autism NMF Project
% Hongyao Ma
% 12/11/2014

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));
load('find_anchor_data_5');

Q_bar = x_colsotc(:, ind_in)';
candidates = 1:size(Q_bar, 1);

%% find_anchors

n_anchor = 10;
[anchor_indices, anchors] = find_anchors(Q_bar, candidates, n_anchor, 1);

figure;
hold on;
scatter3(B(1,:), B(2,:), B(3,:), 'ks');
scatter3(basis(1,:), basis(2,:), basis(3,:), 'ko');
scatter3(Q_bar(:, 1), Q_bar(:, 2), Q_bar(:, 3), 'b.');
scatter3(anchors(:, 1), anchors(:, 2), anchors(:, 3), 'ro');

%% find_anchors_conv

[anchor_indices, anchors] = find_anchors_conv(Q_bar, candidates, n_anchor, 1);

figure;
hold on;
scatter3(B(1,:), B(2,:), B(3,:), 'ks');
scatter3(basis(1,:), basis(2,:), basis(3,:), 'ko');
scatter3(Q_bar(:, 1), Q_bar(:, 2), Q_bar(:, 3), 'b.');
scatter3(anchors(:, 1), anchors(:, 2), anchors(:, 3), 'ro');

