% Demo of the function gen_Doc

% The Autism NMF Project
% Hongyao Ma
% 01/10/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

% The Dirichlet Parameters
alpha = [0.2, 0.3, 0.5];
R = gen_matrix_R(alpha)

A = [1 1 1 0 0 0 0 1 1; ...
       0 1 1 1 1 1 0 0 0; ...
       0 0 0 0 1 1 1 1 1]';
   
% Sizes
V = size(A, 1);
k = size(A, 2);

% Normalize 
A = A ./ repmat(sum(A), V, 1);
Q = A * R * A'
sum(sum(Q))

% Topics
n_Doc = 1000;
topics = drchrnd(alpha, n_Doc);
R_empirical = topics' * topics / n_Doc

% Documents with fixed length
l_Doc = 5;
x = gen_Docs(topics, A, l_Doc, 'fixed');
sum(x, 2)

% Documents with random length
l_Doc = 5;
x = gen_Docs(topics, A, l_Doc, 'random');
sum(x, 2)
mean(sum(x,2))