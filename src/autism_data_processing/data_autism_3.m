% Analyzing Autism Data III --- Cooccurrence Matrix w / w.o. normalization
% The Autism NMF Project

% Hongyao Ma
% Created:  10/31/2014
% Modified: 10/31/2014

clear;
clc;
close all;
% load('/Users/hma/asd_data_nz.mat')
load('E:\My Files\Documents\data/asd_data_nz.mat');

% The row sums
rowsum = sum(data_nz, 2);
figure;
hist(rowsum(rowsum > 1000));

%% Singular Values

[m, n] = size(data_nz);

% 2nd Moment
M2_raw = data_nz'*data_nz;

% Co-occurance matrix
M2_Co_occur = zeros(n);
for i_row = 1:m;
    if mod(i_row, 10) == 0
        i_row
    end
    Q_i = data_nz(i_row, :)'*data_nz(i_row, :);
    Q_i = Q_i / sum(sum(Q_i));
    M2_Co_occur = M2_Co_occur + Q_i;
end

% M2 = M2_raw;
M2 = M2_Co_occur;
figure;
imshow(M2);

%% Normalize row sum

data_nz_rowStoc = data_nz ./ repmat(rowsum, 1, n);
M2_rowStoc = data_nz_rowStoc'*data_nz_rowStoc;

figure;
imshow(M2_rowStoc);


%%  
diagM_TM = zeros(n, 1);

data_modified = zeros(size(data_nz));



%%


% SVD --- saved in the mat and do not have to run
[U, S, V] = svd(M2);

% Plot the singular values
sigmas = diag(S);
r = sum(sigmas > 1e-4)
figure;
semilogy(1:length(sigmas), sigmas+1);
title('log (sigma + 1)');

% Low-rank approximation
i_rank = 50;
sigma1 = sigmas;
sigma1(i_rank+1:end) = 0;
S1 = diag(sigma1);
rankrapprox = U*S1*V';
approx_error = norm(M2 - rankrapprox, 'fro')/norm(M2, 'fro')


%% Low-rank approximations

Sigma2 = zeros(m, n);
error_array = zeros(r, 1);
rankrapprox = zeros(size(M2));

rank_array = 1:100;

for i_rank = 1:length(rank_array)
    if mod(i_rank, 10) == 0
        i_rank
    end
    rankrapprox = rankrapprox + U(:, rank_array(i_rank))*sigmas(rank_array(i_rank)) * V(:, rank_array(i_rank))';
    error_array(rank_array(i_rank)) = norm(M2 - rankrapprox, 'fro')/norm(M2, 'fro');
end

% Plot the errors of low-rank approximations
figure;
plot(1:r, error_array);

stem(1:100, error_array(1:100));


