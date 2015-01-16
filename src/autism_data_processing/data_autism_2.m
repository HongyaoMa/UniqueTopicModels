% Analyzing Autism Data II --- Singular Values and Low-rank Approximation
% The Autism NMF Project

% Hongyao Ma
% Created:  10/20/2014
% Modified: 10/20/2014

clear;
clc;
close all;
% load('/Users/hma/asd_data_nz.mat')
load('E:\My Files\Documents\data/asd_data_nz.mat');

%% Singular Values

[n_patient_nz, n_code_nz] = size(data_nz);

% SVD --- saved in the mat and do not have to run
% [U, S, V] = svd(data_nz);

% Plot the singular values
sigmas = diag(S);
r = sum(sigmas > 1e-4)
figure;
semilogy(1:length(sigmas), sigmas+1);
title('log (sigma + 1)');

% Low-rank approximation
r_low = 50;
sigma1 = sigmas;
sigma1(r_low+1:end) = 0;
S1 = [diag(sigma1); zeros(n_patient_nz - n_code_nz, n_code_nz)];
rankrapprox = U*S1*V';
approx_error = norm(data_nz - rankrapprox, 'fro')/norm(data_nz, 'fro')


%% Low-rank approximations

Sigma2 = zeros(n_patient_nz, n_code_nz);
error_array = zeros(r, 1);
rankrapprox = zeros(size(data_nz));

for r_low = 1:r
    if mod(r_low, 50) == 0
        r_low
    end
    rankrapprox = rankrapprox + U(:,r_low) * sigmas(r_low) * V(:, r_low)';
    error_array(r_low) = norm(data_nz - rankrapprox, 'fro')/norm(data_nz, 'fro');
end

% Plot the 
figure;
plot(1:r, error_array);

%% figure;
stem(1:100, error_array(1:100));

