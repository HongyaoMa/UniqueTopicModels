% Demo of the function compute_err

% The Autism NMF Project
% Hongyao Ma
% 01/16/2015

% clear;
% clc;
% close all;


%% l1 and l2 errors
x1 = [1 2 3; 4 5 6]
x2 = [1 2 2; 4 7 6]

l1_err = compute_err(x1, x2, 'l1')
l2_err = compute_err(x1, x2, 'l2')

% Default
l2_by_default = compute_err(x1, x2)

%% KL Divergence
p1 = [0.2, 0.3, 0.5; 0.1, 0.2, 0.7]
p2 = [0.25, 0.25 0.5; 0.1, 0.2, 0.7]

l2_KL = compute_err(p1, p2, 'row_KL')

% %% Errors
% 
% % Undefined with zero elements in p1
% p1 = [0, 0.5, 0.5; 0.1, 0.2, 0.7]
% p2 = [0.25, 0.25 0.5; 0.1, 0.2, 0.7]
% 
% l2_KL = compute_err(p1, p2, 'row_KL')
% 
% % Do not sum up to 1
% p1 = [0, 0.5, 0.5; 0.1, 0.2, 0.7]
% p2 = [0.25, 0.25 0.5; 0.1, 0.2, 0.7]
% 
% l2_KL = compute_err(p1, p2, 'row_KL')
