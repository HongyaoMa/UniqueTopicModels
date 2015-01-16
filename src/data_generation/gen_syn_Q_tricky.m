% Topic distributions that are difficult to distinguish
% The Autism NMF Project

% Hongyao Ma
% 12/31/2014

clear;
clc;
close all;

%% Topic Distributions

% Anchors: 1, 2, 3?
A_1 = repmat(eye(3), 1, 3);
% A_1 = [1 0 0 1 0 0 1 0 0; ...
%         0 1 0 0 1 0 0 1 0; ...
%         0 0 1 0 0 1 0 0 1];

% Anchors: 1, 4, 7
% Correctly recovered both anchors and topic distributions
A_2 = [1 1 1 0 0 0 0 1 1; ...
       0 1 1 1 1 1 0 0 0; ...
       0 0 0 0 1 1 1 1 1];

% Anchors: 1, 4, 7
% Correctly recovered both anchors and topic distributions
A_3 = [2 1 1 0 1 1 0 0 0; ...
       0 0 0 2 1 1 0 1 1; ...
       0 1 1 0 0 0 2 1 1 ];

% No clear anchors
A_4 = [1 1 1 1 1 1 0 0 0; ...
       1 1 1 0 0 0 1 1 1; ...
       0 0 0 1 1 1 1 1 1];
   
A_5 = [eye(3), ones(3, 6)];

A = A_5';

% Sizes
V = size(A, 1);
k = size(A, 2);

% Normalize 
A = A ./ repmat(sum(A), V, 1);

%% Dirichlet Parameters
alpha0 = 1;

alpha = [0.31, 0.33, 0.36];
alpha = ones(1, k);

alpha = alpha / sum(alpha) * alpha0;

%% The W Matrix
W = alpha' * alpha ;
W = W + diag(alpha);
W = W/(alpha0 * (alpha0 + 1));

%% The Q Matrix
Q = A * W * A'

save('syn_Q_tricky5');
