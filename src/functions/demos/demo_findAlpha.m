% Demo of find_Alpha algorithm

% The Autism-NMF Project
% Hongyao Ma
% 01/05/2015

clear;
clc;
close all;

alpha = [0.2 0.5 0.3]

alpha0 = 2;

R = gen_MatrixR(alpha, alpha0)

alpha_rec = find_Alpha(R)
