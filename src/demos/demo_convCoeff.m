% Demo of convSpan and inConvHull

% The Autism-NMF Project
% Hongyao Ma
% 12/11/2014

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));

%% Initialization

% Generate Random Points within [0, 1]^p
p = 3;
N = 2000;
x = rand(p, N);
scatter3D(x);

% Project onto the unit simplex
x_colsotc = colStoc(x);
scatter3D(x_colsotc);
rank_xcolStoc = rank(x_colsotc);

%% Basis

B = eye(p);
B1 = [0.8, 0.1, 0.1; 0.1, 0.8, 0.1; 0.1, 0.1, 0.8];
B2 = [0.9, 0.05, 0.05; 0.05, 0.9, 0.05; 0.05, 0.05, 0.9];
B3 = [[0.7; 0.2; 0.1], [0.7; 0.1; 0.2], [0.1; 0.7; 0.2], [0.2; 0.7; 0.1], [0.1; 0.2; 0.7], [0.2; 0.1; 0.7]];
B4 = [[2/3; 1/3; 0], [2/3; 0; 1/3], [0; 2/3; 1/3], [1/3; 2/3; 0], [0; 1/3; 2/3], [1/3; 0; 2/3]];
B5 = [[1/2; 2.5/8; 1.5/8], [1/2; 1.5/8; 2.5/8], [1.5/8; 1/2; 2.5/8], [2.5/8; 1/2; 1.5/8], [1.5/8; 2.5/8; 1/2], [2.5/8; 1.5/8; 1/2]];

% Plot the basis
figure;
hold on;
scatter3(B(1,:), B(2,:), B(3,:), 'k*');
scatter3(B1(1,:), B1(2,:), B1(3,:), 'bo');
scatter3(B2(1,:), B2(2,:), B2(3,:), 'ro');
scatter3(B3(1,:), B3(2,:), B3(3,:), 'go');
scatter3(B4(1,:), B4(2,:), B4(3,:), 'mo');
scatter3(B5(1,:), B5(2,:), B5(3,:), 'ko');

%% Demo of convSpan

% Take the basis B3
basis = B2;
[ind, num] = convSpan(x_colsotc, basis);

figure;
hold on;
scatter3(basis(1,:), basis(2,:), basis(3,:), 'k*');
scatter3(x_colsotc(1,ind), x_colsotc(2,ind), x_colsotc(3,ind), 'b.');
ind_out = logical(1-ind);
scatter3(x_colsotc(1,ind_out), x_colsotc(2,ind_out), x_colsotc(3,ind_out), 'r.');

%% Demo of inConvHull
[ind_in, res_norm] = inConvHull(basis, x_colsotc);

figure;
hold on;
scatter3(basis(1,:), basis(2,:), basis(3,:), 'k*');
scatter3(x_colsotc(1, ind_in), x_colsotc(2, ind_in), x_colsotc(3, ind_in), 'b.');
ind_out = logical(1 - ind_in);
scatter3(x_colsotc(1, ind_out), x_colsotc(2, ind_out), x_colsotc(3, ind_out), 'r.');

figure;
hist(res_norm);

Q_bar = x_colsotc(:, ind_in)';
candidates = 1:size(Q_bar, 1);

save('find_anchor_data_2');




