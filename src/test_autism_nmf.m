% Testing the NMF algorithms on Autism Data
% The Autism NMF Project

% Hongyao Ma
% Created:  10/20/2014
% Modified: 10/20/2014

clear;
clc;
close all;
load('/Users/hma/asd_data_nz.mat')

% Add all subfolders of the source directory
addpath(genpath('../src'));

%% NMF
r = 3;

maxiter = 100000;
timelimit = 60;
% data = data_nz;

non_zero_in_col = sum(data_nz > 0);
data_nz = data_nz(:, non_zero_in_col >= 4);

[U, S, V] = svd(data_nz);

% data_zeropattern = double(data_nz > 0);
% [U, S, V] = svd(data_zeropattern);
% data = data_zeropattern;

U0 = U(:, 1:r);
V0 = V(:, 1:r)';

[Uha,Vha,eha,tha] = HALSacc(data_nz, U0,V0,0.5,0.1,maxiter,timelimit);
disp(sprintf('accelerated HALS terminated with final error %f',eha(end)));
figure;
plot(tha, eha)


%%

data_rec = Uha * Vha;
figure;
imshow(data_rec)


fnorm = norm(data_nz - Uha * Vha, 'fro')
fnorm_normalized = fnorm/norm(data_nz, 'fro')


%% Normalization

Uha1 = bsxfun( @rdivide , Uha , sum( Uha ) );
Vha1 = bsxfun( @rdivide , Vha , sum( Vha ) );

figure;
plot(Uha1);

figure;
plot(Vha1')
