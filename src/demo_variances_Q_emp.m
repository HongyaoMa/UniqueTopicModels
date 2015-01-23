% Compute the variances of the Empirical Q Matrices

% The Autism NMF Project
% Hongyao Ma
% Created:  01/19/2015
% Modified: 01/22/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

%% Parameter Generation

[A, V, k] = gen_trickyA(5);
alpha = gen_alpha(0.3, k, 'uniform');

% % Set 2
% V = 500;
% k = 5;
% p_anchor = 0.05;
% A = gen_matrix_A(V, k, p_anchor);
% alpha = gen_alpha(0.03, k, 'uniform');

R = gen_matrix_R(alpha);
Q = A * R * A';
Q_bar = rowStoc(Q);

%% Single Documents

n_Doc = 1000;
l_Doc = 20;

type_normalize = 'none';
type_length = 'fixed';

topics = drchrnd(alpha, n_Doc);
x = gen_Docs(topics, A, l_Doc, type_length);
        
var1 = zeros(V, V);
var2 = zeros(V, V);
var = zeros(V,V);

for i_Doc = 1:n_Doc
    M2 = (x(i_Doc, :)' * x(i_Doc, :) - diag(x(i_Doc, :))) / l_Doc / (l_Doc - 1);    
    M2_exp = A * topics(i_Doc, :)' * topics(i_Doc, :) * A';
    var1 = var1 + (M2 - M2_exp).^2;
    var2 = var2 + (Q - M2_exp).^2;  
    var = var + (M2 - Q).^2;    
end

%% Results
% The Mean's
sumvar = norm(var, 'fro') 
sumvar1 = norm(var1, 'fro') 
sumvar2 = norm(var2, 'fro') 

err = norm(var1 + var2 - var, 'fro')/ norm(var, 'fro')

% Visualization
figure;
subplot(2,2,1)
title('Q')
hold on;
imagesc(Q);

subplot(2,2,2)
title('E[(xx^T - diag x - Q)^2]')
hold on;
imagesc(var)

subplot(2,2,3)
title('E[(xx^T - diag x - Aww^TA^T)^2]')
hold on;
imagesc(var1)

subplot(2,2,4)
title('E[(Aww^TA^T - Q)^2]')
hold on;
imagesc(var2)

colormap('gray')


%% Sets of documents
n_iter = 5000;
n_Doc = 10000;
l_Doc = 20;

type_length = 'fixed';

var = zeros(V);
var1 = zeros(V);
var2 = zeros(V);
tic
for i_iter = 1:n_iter
    if mod(i_iter, 5000) == 0
        i_iter
        toc
    end
    topics = drchrnd(alpha, n_Doc);
    x = gen_Docs(topics, A, l_Doc, type_length);
    [Q_emp, Q_emp_bar] = gen_matrix_Q(x, 0, 'none');
    
    R_emp = topics' * topics / n_Doc;    
    M2_exp = A * R_emp * A';

    var = var + (Q_emp - Q).^2;   
    var1 = var1 + (Q_emp - M2_exp).^2;
    var2 = var2 + (Q - M2_exp).^2;  
end

%% The Mean's
sumvar = norm(var, 'fro')
sumvar1 = norm(var1, 'fro')
sumvar2 = norm(var2, 'fro') 

err = norm(var1 + var2 - var, 'fro')/ norm(var, 'fro')

% Visualization
figure;
subplot(2,2,1)
title('Q')
hold on;
imagesc(Q);

subplot(2,2,2)
title('E[(Q_{emp} - Q)^2]')
hold on;
imagesc(var)

subplot(2,2,3)
title('E[(Q_{emp} - A*R_{emp}*A^T)^2]')
hold on;
imagesc(var1)

subplot(2,2,4)
title('E[( A*R_{emp}*A^T - Q)^2]')
hold on;
imagesc(var2)

colormap('gray')

%% Analytical result
p = sum(Q)';
var_analytical = var_Q(p, l_Doc);

% For a set of n_Doc documents:
var_analytical = var_analytical / n_Doc;

err_var = var1/n_iter - var_analytical;

err = norm(err_var, 'fro') / norm(var1/n_iter, 'fro')
err_diag = norm(diag(err_var)) / norm(diag(var1/n_iter))

figure;
subplot(1,2,1);
imagesc(var1);

subplot(1,2,2)
imagesc(err_var);

colormap('gray')

% No Diagonal
err_var_nodiag = err_var;
var1_nodiag = var1/n_iter;
var_analytical_nodiag = var_analytical;
 
for i = 1:V
    err_var_nodiag(i,i) = 0;
    var1_nodiag(i,i) = 0;
    var_analytical_nodiag(i,i) = 0;
end

err_nodiag = norm(err_var_nodiag, 'fro') / norm(var1_nodiag, 'fro')

figure;
subplot(1,2,1);
imagesc(var1_nodiag);

subplot(1,2,2)
imagesc(var_analytical_nodiag);

colormap('gray')

% for i = 1:V
%     var_analytical(i,i) = 0;
%     for j = 1:V
%         var_analytical(i,j) = ...
%             (-4 * p(i)^2 * p(j)^2 + p(i)^2*p(j) + p(j)^2 *p(i)) * l_Doc^3 + ...
%             (10*p(i)^2* p(j)^2 - 3*(p(i)^2*p(j) + p(j)^2 *p(i)) + p(i) *p(j)) * l_Doc^2 + ...
%             (6*p(i)^2* p(j)^2 + p(i)^2*p(j) + p(j)^2 *p(i) - p(i) *p(j))*l_Doc;
%     end
% end

%%
err_var = err_var + var2/n_iter;
err = norm(err_var, 'fro') / norm(var1/n_iter, 'fro')

