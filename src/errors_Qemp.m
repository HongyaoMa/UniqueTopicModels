% Demo of l1-error in Q_emp

% The Autism NMF Project
% Hongyao Ma
% Created:  01/11/2015
% Modified: 01/16/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

%% Distributions

% Set 1
[A, V, k] = gen_trickyA(2);
alpha = gen_alpha(0.3, k, 'uniform');

% % Set 2
% V = 500;
% k = 5;
% p_anchor = 0.05;
% A = gen_matrix_A(V, k, p_anchor);
% alpha = gen_alpha(0.3, k, 'random');

R = gen_matrix_R(alpha);
% R = [1 0 0; 0 0.45 0.55; 0 0.55 0.45];

Q = A * R * A';
Q_bar = rowStoc(Q);

%%

N_iter = 2000;
% Topics
n_Doc = 500;
l_Doc = 20;

flag_no_normalize = 1;

err_im_l1 = zeros(size(Q));
err_im_bar_l1 = zeros(size(Q));

err_im_l2 = zeros(size(Q));
err_im_bar_l2 = zeros(size(Q));

err_Q = zeros(1, N_iter);
err_Q_anchor = zeros(1, N_iter);

err_im_rand = zeros(size(Q));
err_im_bar_rand = zeros(size(Q));

err_Q_rand = zeros(1, N_iter);
err_Q_anchor_rand = zeros(1, N_iter);

err_Q_rand_none = zeros(1, N_iter);
err_Q_rand_n32 = zeros(1, N_iter);
err_Q_rand_n = zeros(1, N_iter);

tic
for i = 1:N_iter
    
    % Where are we?
    if mod(i, 100) == 0
        i
        toc
    end
    
    topics = drchrnd(alpha, n_Doc);
    R_empirical = topics' * topics / n_Doc;
    
    %% Documents with fixed length
    x = gen_Docs(topics, A, l_Doc, 'fixed');
    
    % Word-Word Co-occurrance
    [Q_emp, Q_emp_bar] = gen_matrix_Q(x, 0, 'original');
    
    % Normalized Squared Error
    err_Q(i) = norm(Q - Q_emp, 'fro') / norm(Q, 'fro');
    err_Q_anchor(i) = norm(Q(1:3,:) - Q_emp(1:3,:), 'fro') / norm(Q(1:3,:), 'fro');
    
    err_im_l1 = err_im_l1 + abs(Q - Q_emp);
    err_im_bar_l1 = err_im_bar_l1 + abs(Q_bar - Q_emp_bar);
    
    err_im_l2 = err_im_l2 + (Q - Q_emp).^2;
    err_im_bar_l2 = err_im_bar_l2 + (Q_bar - Q_emp_bar).^2;
    
    %% Documents with random length
    x_rand = gen_Docs(topics, A, l_Doc, 'random');
    
    % Word-Word Co-occurrance
    [Q_emp_rand, Q_emp_bar_rand] = gen_matrix_Q(x_rand, 0, 'original');
    
    % Normalized Squared Error
    err_Q_rand(i) = norm(Q - Q_emp_rand, 'fro') / norm(Q, 'fro');
    err_Q_anchor_rand(i) = norm(Q(1:3,:) - Q_emp_rand(1:3,:), 'fro') / norm(Q(1:3,:), 'fro');
    
    err_im_rand = err_im_rand + abs(Q - Q_emp_rand);
    err_im_bar_rand = err_im_bar_rand + abs(Q_bar - Q_emp_bar_rand);
    
    % No Normalization
    [Q_emp_rand_none, ~] = gen_matrix_Q(x_rand, 0, 'none');
    err_Q_rand_none(i) = norm(Q - Q_emp_rand_none, 'fro') / norm(Q, 'fro');    
    
    [Q_emp_rand_n32, ~] = gen_matrix_Q(x_rand, 0, 'n32');
    err_Q_rand_n32(i) = norm(Q - Q_emp_rand_n32, 'fro') / norm(Q, 'fro');
    
    [Q_emp_rand_n, ~] = gen_matrix_Q(x_rand, 0, 'n32');
    err_Q_rand_n(i) = norm(Q - Q_emp_rand_n, 'fro') / norm(Q, 'fro');  
    
end

% Normalize the errors
err_im_l1 = err_im_l1 / N_iter;
err_im_bar_l1 = err_im_bar_l1 / N_iter;

err_im_l2 = err_im_l1 / N_iter;
err_im_bar_l2 = err_im_bar_l1 / N_iter;

%% Mean errors
mean_fixedlength = mean(err_Q)
mean(err_Q_anchor);

mean_original = mean(err_Q_rand)
mean(err_Q_anchor_rand);

mean_none = mean(err_Q_rand_none)
mean_n32 = mean(err_Q_rand_n32)
mean_n = mean(err_Q_rand_n)

%% Visualize the l1 errors

n_imshow = 5;
figure;

subplot(2, 2, 1);
hold on;
title('average(|Q - Q_{emp}|)')
imshow(err_im_l1 * 50 * n_imshow);

subplot(2, 2, 2);
hold on;
title('average(|Q - Q_{emp}|) ./ Q')
imshow(err_im_l1 ./ Q  * n_imshow);

subplot(2, 2, 3);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|)')
imshow(err_im_bar_l1 *10 * n_imshow);

subplot(2, 2, 4);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|) ./ Q_{bar}')
imshow(err_im_bar_l1./ Q_bar  * n_imshow);

colormap('gray')



%% Compare the fixed lengths and the random lengths
figure;

subplot(2, 2, 1);
hold on;
title('average(|Q - Q_{emp}|)')
imshow(err_im_l1 * 50 * n_imshow);

subplot(2, 2, 2);
hold on;
title('average(|Q - Q_{emp, rand}|)')
imshow(err_im_rand * 50 * n_imshow);

subplot(2, 2, 3);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|)')
imshow(err_im_bar_l1 *10 * n_imshow);

subplot(2, 2, 4);
hold on;
title('average(|Q_{bar} - Q_{bar, emp, rand}|)')
imshow(err_im_bar_rand *10 * n_imshow);


%% Visualize the l2 errors
n_imshow = 10000;
figure;

subplot(2, 2, 1);
hold on;
title('average(|Q - Q_{emp}|)')
imshow(err_im_l2 * 50 * n_imshow);

subplot(2, 2, 2);
hold on;
title('average(|Q - Q_{emp}|) ./ Q')
imshow(err_im_l2 ./ Q  * n_imshow);

subplot(2, 2, 3);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|)')
imshow(err_im_bar_l2 *10 * n_imshow);

subplot(2, 2, 4);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|) ./ Q_{bar}')
imshow(err_im_bar_l2./ Q_bar  * n_imshow);

colormap('gray')