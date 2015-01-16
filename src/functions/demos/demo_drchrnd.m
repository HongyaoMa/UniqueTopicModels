% Demo of the function drchrnd

% The Autism NMF Project
% Hongyao Ma
% 01/06/2015

clear;
clc;
close all;

addpath(genpath('../src'));

% Parameters
alpha = [2 3 5];

% Check the summation = 1
n = 5;
data = drchrnd(alpha, n)
sum(data, 2)

% Mean = alpha_i / alpha_0
n = 1000;
data = drchrnd(alpha, n)
mean(data)

% Covariance Matrix
R = gen_MatrixR(alpha)
R_empirical = data' * data / n;

