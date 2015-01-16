% Analyzing the CBCL Dataset
% The NMF Autism Project

% Hongyao Ma
% 10/21/2014

clear;
clc;
close all;

load cbclim; 
[U, S, V] = svd(M);

%% Singular Values

[m, n] = size(M);

% Plot the singular values
sigmas = diag(S);
r = sum(sigmas > 1e-4)
figure;
semilogy(1:length(sigmas), sigmas+1);
title('log (sigma + 1)');


%% Low-rank approximations

Sigma2 = zeros(m, n);
error_array = zeros(r, 1);
rankrapprox = zeros(size(M));

for r_low = 1:r
    if mod(r_low, 50) == 0
        r_low
    end
    rankrapprox = rankrapprox + U(:,r_low)*sigmas(r_low) * V(:, r_low)';
    error_array(r_low) = norm(M - rankrapprox, 'fro')/norm(M, 'fro');
end

% Plot the 
figure;
plot(1:r, error_array);

%% figure;
stem(1:100, error_array(1:100));

%% Heatmap
figure;
imshow(M);