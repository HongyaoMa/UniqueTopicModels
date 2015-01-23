% Analyzing the Q_matrices
% Find how many rows are covered by those that correspond to words with
% certain dpw
% The method that determines "covered" is not correct. However this doesn't
% really matter since the Q matrix is generally full rank
% The other method using convex span is too slow to run.

% The NMF Project

% Hongyao Ma
% Created:  11/25/2014
% Updated:  12/11/2014

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));

load('data_Q_10_10.mat');

%% 
% How many docs in which a word showed up
dpw = sum(data_new > 0);

% Testing for an arrway of niminum # of docs
max_dpw = 100;
dpw_array = 0:5:max_dpw;
n_dpw = length(dpw_array);

% Initialization of counters
num_covered = zeros(n_dpw, 1);
num_big_dpw = zeros(n_dpw, 1);
mean_err = zeros(n_dpw, 1);

for i_dpw = 1:n_dpw
    
    % indices of words with dpw at least dpw_array(i_dpw)
    ind_big_dpw = dpw >= dpw_array(i_dpw);
    
    % Number of words showing up in at least dpw_array(i_dpw) docs
    num_big_dpw(i_dpw) = sum(ind_big_dpw);
    
    % Extract corresponding rows in the conditional Q matrix
    Q_dpw = Q_bar(ind_big_dpw, :);
    
    % Updated Method --- too slow to run
    [ind_in, res_norm] = inConvHull(Q_dpw', Q_bar');
    num_covered(i_dpw) = sum(ind_in);
    mean_err(i_dpw) = mean(res_norm);
    
%     % Old Method --- not really accurate but should be fine
%     % Compute the spanning coefficients
%     coeff = (Q_dpw') \ Q_bar';
%     
%     % number of rows with non-negative coefficients
%     num_covered(i_dpw) = sum(all(coeff > -1e-8));   
%     
%     % The non-negative coefficients
%     coeff_nonneg = coeff;
%     coeff_nonneg(coeff < 0) = 0;
%     
%     % TODO
%     Q_bar_hat = ((Q_dpw') * coeff_nonneg)';
%     error = sum((Q_bar_hat - Q_bar).^2, 2) ./ sum(Q_bar.^2, 2);
%     size(error)
%     mean_err(i_dpw) = mean(error);

end

figure;
stem(dpw_array, num_big_dpw, 'k*');
hold on;
stem(dpw_array, num_covered);

figure;
stem(dpw_array, mean_err);
