% Demo of the lsqlin function

% The Autism NMF Project
% Hongyao Ma
% 12/11/2014

clear;
clc;
close all;

B1 = [0.8, 0.1, 0.1; 0.1, 0.8, 0.1; 0.1, 0.1, 0.8];
C = B1;
d = [0.9, 0.05, 0.05]';
A = [];
b = [];

Aeq = [1, 1, 1];
beq = 1;
lb = [0, 0, 0]';
ub = [1, 1, 1]';

[x, resnorm, residual, exitflag, output, lambda] = lsqlin(C, d, A, b, Aeq, beq, lb, ub)
