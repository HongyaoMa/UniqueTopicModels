% Analyzing Autism Data V --- Cleaning the data, construct the Q matrix,
% save for later use and do the low-rank approximation
% The Autism NMF Project

% Hongyao Ma
% Created:  11/24/2014
% Modified: 11/24/2014

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));

%% Data
% Load the data with no all zero rows/columns
% load('/Users/hma/asd_data_nz.mat')
load('E:\My Files\Documents\data/asd_data_nz.mat');

% Clean the data
data_new = clean_data(data_nz, [10 0 10 0]);

%% Co-occurrence Matrices
[Q, Q_bar] = gen_matrix_Q(data_new, 1);

% Save the data to be loaded by Python
M = sparse(data_new');
save('data_10_10', 'M');

% Save the data for later use
save('data_Q_10_10', 'data_new', 'Q', 'Q_bar');

%% Low Rank Approximation

% SVD and Singular Values
[U S V] = svd(Q);

figure;
semilogy(1:length(diag(S)), diag(S), 'o');
axis([0, length(diag(S)), 1e-10, 1])
xlabel('# of Singular Value');

% Rank 50 Approximation
S50 = S;
for i = 51:length(diag(S));
    S50(i,i) = 0;
end;

Q50 = U * S50 * V';
error = norm(Q50 - Q, 'fro')/norm(Q, 'fro')

%% Show the figure and the difference

figure;
imshow(Q50 * size(data_new,1) ^2);
title('Q50');

Q_diff = Q - Q50;
figure;
imshow(Q_diff * size(data_new,1)^2);
title('Q - Q50');

%% Low-Rank Approximation

r = 1:100;
err_r = zeros(size(r));

Q_r = zeros(size(Q));

for i = 1:length(r)
    Q_r = Q_r + U(:, r(i))* S(r(i), r(i)) * V(:, r(i))';
    err_r(i) = norm(Q_r - Q, 'fro')/norm(Q, 'fro');
end

figure;
stem(err_r);
xlabel('rank');
ylabel('||Q - Q_r||_F / ||Q||_F');



