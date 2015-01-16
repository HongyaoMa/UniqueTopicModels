% Compute the distances from convex hull
% The NMF Autism Project

% Hongyao Ma
% 01/05/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));

B = [0.8, 0.1, 0.1; ...
     0.1, 0.8, 0.1; ...
     0.1, 0.1, 0.8]

[ind_in, res_norm, coeffs] = inConvHull(B(:, 1:2), B(:,3))

b1 = B(:,1);
b2 = B(:,2);
b3 = B(:,3);

c1 = 0.2;
c2 = 0.4;
c3 = 0.4;

x = b1*c1 + b2 * c2 + b3*c3;

[ind_in, res_norm, coeffs] = inConvHull(B(:, 1:2), x)
