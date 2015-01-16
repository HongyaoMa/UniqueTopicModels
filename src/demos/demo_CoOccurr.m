% Demo Outer Product
% The Autism Project

% Hongyao Ma
% 11/07/2014

clear;
clc;
close all;

p1 = [1 3 5];
p2 = [2 2 2];

X = [p1; p2];

X = X./ repmat(sum(X, 2), 1, 2);

X'*X

p1'*p1 + p2'*p2
