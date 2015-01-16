% Demo of moments of Multinomial Distributed Variables

% The Autism NMF Project
% Hongyao Ma
% 01/15/2015

clear;
clc;
close all;

% The probability
p = [0.2 0.3 0.5];
R = p'*p

% # of instances
N_iter = 1000;

n1 = 1;
x1 = mnrnd(n1, p, N_iter);
x1_2 = mnrnd(n1, p, N_iter);
x1_mean = mean(x1);

R1 = x1' * x1_2 / N_iter

n2 = 2;
x2 = mnrnd(n2, p, N_iter);
x2_mean = mean(x2)/n2;

R2 = (x2' * x2 - diag(sum(x2)))/ N_iter