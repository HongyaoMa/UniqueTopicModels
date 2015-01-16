% Demo of plot_SingularVals, err_low_approx, and err_rank_r_approx

% The Autism NMF Project
% Hongyao Ma
% 01/11/2015

clear;
clc;
close all;

%% Generate the data matrix
% The Dirichlet Parameters
alpha = [0.2, 0.3, 0.5];
R = gen_matrix_R(alpha);

A = [1 1 1 0 0 0 0 1 1; ...
       0 1 1 1 1 1 0 0 0; ...
       0 0 0 0 1 1 1 1 1]';
   
% Sizes
V = size(A, 1);
k = size(A, 2);

% Normalize 
A = A ./ repmat(sum(A), V, 1);
Q = A * R * A';

% Topics
n_Doc = 2000;
topics = drchrnd(alpha, n_Doc);
R_empirical = topics' * topics / n_Doc;

% Documents
l_Doc = 50;
x = gen_Docs(topics, A, l_Doc);

% Word-Word Co-occurrance
[Q_emp, Q_bar_emp] = gen_matrix_Q(x, 0);

%% Demo of the plot_SingularVals
plot_SingularVals(Q)
plot_SingularVals(Q_emp)

%% Demo of the err_low_r_approx
err_array = low_rank_approx(Q_emp, 1:length(Q), 1)
err = rank_r_approx(Q_emp, 5)
