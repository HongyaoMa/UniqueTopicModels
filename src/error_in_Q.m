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
% Anchors: 1, 2, 3?
A_1 = repmat(eye(3), 1, 3);
% A_1 = [1 0 0 1 0 0 1 0 0; ...
%         0 1 0 0 1 0 0 1 0; ...
%         0 0 1 0 0 1 0 0 1];

% Anchors: 1, 4, 7
% Correctly recovered both anchors and topic distributions
A_2 = [1 1 1 0 0 0 0 1 1; ...
    0 1 1 1 1 1 0 0 0; ...
    0 0 0 0 1 1 1 1 1];

% Anchors: 1, 4, 7
A_3 = [2 1 1 0 1 1 0 0 0; ...
    0 0 0 2 1 1 0 1 1; ...
    0 1 1 0 0 0 2 1 1 ];

% No clear anchors
A_4 = [1 1 1 1 1 1 0 0 0; ...
    1 1 1 0 0 0 1 1 1; ...
    0 0 0 1 1 1 1 1 1];

A_5 = [eye(3), ones(3, 6)];

A_6 = [eye(3)*3, eye(3)*2, eye(3)];
A_7 = [eye(3)*9, eye(3)*3, eye(3)];

A = A_1';
A = A_2';

% Sizes
V = size(A, 1);
k = size(A, 2);

% Normalize
A = A ./ repmat(sum(A), V, 1);

alpha = ones(1, k);
R = gen_matrix_R(alpha, 1);
% R = [1 0 0; 0 0.45 0.55; 0 0.55 0.45];

Q = A * R * A';
Q_bar = rowStoc(Q);

%%

N_iter = 2000;
% Topics
n_Doc = 500;
l_Doc = 20;

flag_no_normalize = 1;

err_im = zeros(size(Q));
err_im_bar = zeros(size(Q));

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
    % [Q_emp, Q_emp_bar] = gen_matrix_Q(x, 0)
    [Q_emp, Q_emp_bar] = gen_matrix_Q(x, 0, 'original');
    
    % Normalized Squared Error
    err_Q(i) = norm(Q - Q_emp, 'fro') / norm(Q, 'fro');
    err_Q_anchor(i) = norm(Q(1:3,:) - Q_emp(1:3,:), 'fro') / norm(Q(1:3,:), 'fro');
    
    err_im = err_im + abs(Q - Q_emp);
    err_im_bar = err_im_bar + abs(Q_bar - Q_emp_bar);
    
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


%% Mean errors
mean_fixedlength = mean(err_Q)
mean(err_Q_anchor);

mean_original = mean(err_Q_rand)
mean(err_Q_anchor_rand);

mean_none = mean(err_Q_rand_none)
mean_n32 = mean(err_Q_rand_n32)
mean_n = mean(err_Q_rand_n)

%% Visualize the errors

n_imshow = 1;
figure;

subplot(2, 2, 1);
hold on;
title('average(|Q - Q_{emp}|)')
imshow(err_im/N_iter * 50 * n_imshow);

subplot(2, 2, 2);
hold on;
title('average(|Q - Q_{emp}|) ./ Q')
imshow(err_im ./ Q / N_iter * n_imshow);

subplot(2, 2, 3);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|)')
imshow(err_im_bar / N_iter*10 * n_imshow);

subplot(2, 2, 4);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|) ./ Q_{bar}')
imshow(err_im_bar./ Q_bar / N_iter * n_imshow);

colormap('gray')

%% Compare the fixed lengths and the random lengths
figure;

subplot(2, 2, 1);
hold on;
title('average(|Q - Q_{emp}|)')
imshow(err_im/N_iter * 50 * n_imshow);

subplot(2, 2, 2);
hold on;
title('average(|Q - Q_{emp, rand}|)')
imshow(err_im_rand/N_iter * 50 * n_imshow);

subplot(2, 2, 3);
hold on;
title('average(|Q_{bar} - Q_{bar, emp}|)')
imshow(err_im_bar / N_iter*10 * n_imshow);

subplot(2, 2, 4);
hold on;
title('average(|Q_{bar} - Q_{bar, emp, rand}|)')
imshow(err_im_bar_rand / N_iter*10 * n_imshow);


% imshow([err_im/N_iter*50, err_im ./ Q / N_iter; err_im_bar/N_iter*10, err_im_bar./ Q_bar / N_iter ]*10 );
%
% figure;
% imshow(err_im ./ Q / 2000);

% % Histogram of the errors
% figure;
% hist(err_Q);
% figure;
% hist(err_Q_anchor);
