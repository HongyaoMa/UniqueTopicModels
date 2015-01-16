% Demo of the function gen_matrix_A

% The Autism NMF Project
% Hongyao Ma
% 01/11/2015

clear;
clc;
close all;

% Parameters
V = 9;              % # of words
k = 3;              % # of topics
p_anchor = 0.02;    % Prob. of the anchors

% Generate the matrix A
A = gen_matrix_A(V, k, p_anchor)
sum(A)


% Generate the matrix A with multiple set of anchors
A = gen_matrix_A(V, k, [0.1 0.2])
sum(A)