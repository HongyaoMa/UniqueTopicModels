% Generating Synthetic Q Matrices

% The Autism NMF Project
% Hongyao Ma
% 12/17/2014

clear;
clc;
close all;

% Size of dictionary
V = 500;

% Number of topics
k = 20;

%% The A Matrix
% Probability for Anchors Words
p_anchor = 0.01;
A1 = eye(k) * p_anchor;

A2 = rand(V-k, k);
%A2 = randi(2, V-k, k)-1

sumA2 = sum(A2);
A2 = (1-p_anchor)* A2./ repmat(sumA2, V-k, 1);

A = [A1; A2];
sum(A)

%% Dirichlet Parameters
alpha0 = 0.3;
% alpha = [0.2, 0.3, 0.5];
alpha = rand(1, k);
alpha = alpha / sum(alpha) * alpha0;

%% The W Matrix
R = alpha' * alpha ;
R = R + diag(alpha);
R = R/(alpha0 * (alpha0 + 1));

%% The Q Matrix
Q = A * R * A';

save('syn_Q_500_20_000001');


