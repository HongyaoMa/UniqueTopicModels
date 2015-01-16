% Demo of the Basis-shifting matrices

% The Autism-NMF Project
% Hongyao Ma
% 10/13/2014

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

% Original rank
rankx = rank(x)
[u, s, v] = svds(x);
diag(s)

% Project onto the unit simplex
x_colsotc = colStoc(x);
scatter3D(x_colsotc);
rank_xcolStoc = rank(x_colsotc)

% SVD of the col-stochastic matrix
[ucol, scol, vcol] = svds(x_colsotc);
diag(scol)

%% Basis

% Canonical Basis
B = eye(p);
conv_coeff = rand_simplex(p, N);

x = B * conv_coeff;
scatter3D(x);
hold on;
scatter3(B(:,1), B(:,2), B(:,3), 'k*');

rank_can = rank(x)
[u, s, v] = svds(x);
diag(s)

% Naive Basis 1
B1 = [0.8, 0.1, 0.1; 0.1, 0.8, 0.1; 0.1, 0.1, 0.8];
conv_coeff = rand_simplex(p, N);

x_1 = B1 * conv_coeff;
scatter3D(x_1);
hold on;
scatter3(B1(1,:), B1(2,:), B1(3,:), 'k');
scatter3(B(1,:), B(2,:), B(3,:), 'k*');


% Naive Basis 2
B2 = [0.9, 0.05, 0.05; 0.05, 0.9, 0.05; 0.05, 0.05, 0.9];
conv_coeff = rand_simplex(size(B2, 2), N);

x_2 = B2 * conv_coeff;
scatter3D(x_2);
hold on;
scatter3(B2(1,:), B2(2,:), B2(3,:), 'k');
scatter3(B(:,1), B(:,2), B(:,3), 'k*');

%% Basis Shifting

% Say B1 = B2Q21
Q21 = B2\B1
rankQ21 = rank(Q21)
sumQ21 = sum(Q21)

% B2 = B1Q12
Q12 = B1\B2
rankQ12 = rank(Q12)
sumQ12 = sum(Q12)

% Inverses
Q21inv = Q21^(-1)
Q12inv = Q12^(-1)

%% More generators than the rank
% N = 100000;
B3 = [[0.7; 0.2; 0.1], [0.7; 0.1; 0.2], [0.1; 0.7; 0.2], [0.2; 0.7; 0.1], [0.1; 0.2; 0.7], [0.2; 0.1; 0.7]]

conv_coeff = rand_simplex(size(B3, 2), N);
x_3 = B3 * conv_coeff;

scatter3D(x_3);
hold on;
scatter3(B3(1,:), B3(2,:), B3(3,:), 'k');
scatter3(B(:,1), B(:,2), B(:,3), 'k*');
% TODO: points would concentrate in the middle
% It would be fun to compute the distribution

rankB3 = rank(B3)
rankx3 = rank(x_3)

%% NMF

% This basically gives you B3 = I * B3
[W,H] = nnmf(B3,3)

B4 = [[0.7; 0.2; 0.1], [0.1; 0.7; 0.2], [0.2; 0.1; 0.7]]
coeff = B4\[0.7; 0.1; 0.2]
