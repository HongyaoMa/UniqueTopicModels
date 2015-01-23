% Demo of var1_Q
% The Autism NMF Project
% Hongyao Ma
% 01/22/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

%% Parameter Generation

[A, V, k] = gen_trickyA(5);
alpha = gen_alpha(0.3, k, 'uniform');


R = gen_matrix_R(alpha);
Q = A * R * A';
Q_bar = rowStoc(Q);

%% Single Distribution

n_iter = 1000000;
l_Doc = 20;

type_normalize = 'none';
type_length = 'fixed';

topic = drchrnd(alpha, 1);
topics = repmat(topic, n_iter, 1);
x = gen_Docs(topics, A, l_Doc, type_length);
        
var1 = zeros(V, V);
var2 = zeros(V, V);
var = zeros(V,V);

for i_iter = 1:n_iter
    M2 = (x(i_iter, :)' * x(i_iter, :) - diag(x(i_iter, :))) / l_Doc / (l_Doc - 1);    
    M2_exp = A * topics(i_iter, :)' * topics(i_iter, :) * A';
    var1 = var1 + (M2 - M2_exp).^2;
    var2 = var2 + (Q - M2_exp).^2;  
    var = var + (M2 - Q).^2;    
end

var_analytical = var_Q(A * topic', l_Doc);

err_var = var1/n_iter - var_analytical;

err = norm(err_var, 'fro') / norm(var1/n_iter, 'fro') 
err_diag = norm(diag(err_var)) / norm(diag(var1/n_iter))