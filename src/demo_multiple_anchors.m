% Demo --- Performance of find_anchors and Recover while we have mltiple
% sets of anchors with different probabilities

% The Autism NMF Project
% Hongyao Ma
% Created:  01/12/2015
% Modified: 01/12/2015

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));
addpath(genpath('../data'));

%% Parameters

V = 1000;
k = 5;
n_Doc = 2000;
l_Doc = 500;
n_iter = 10;

% Test - doing no normalization
flag_no_norm = 0;

% Test - adjust the anchors
flag_adjust = 1;

% Test - using true anchors
flag_trueanchors = 1;

% Plot the histograms
flag_histogram = 0;

%% Anchors and Topic Distributions

% Probability of sets of anchors
p_anchor = [0.001, 0.01, 0.1];
n_anchor = length(p_anchor);        % # of sets of anchors
p_anchor = sort(p_anchor);          % Sort w.r.t. probabilities
worst_anchors = 1:k;                % Anchors with lowest probability
best_anchors = (n_anchor-1)*k + 1 : n_anchor*k; % with highest probability

% Topic distribution matrix A
A = gen_matrix_A(V, k, p_anchor);

% %%%%%%%
% A = [1 1 1 0 0 0 0 1 1; ...
%      0 1 1 1 1 1 0 0 0; ...
%      0 0 0 0 1 1 1 1 1]';
%
% % Sizes
% V = size(A, 1);
% k = size(A, 2);
%
% % Normalize
% A = A ./ repmat(sum(A), V, 1);
% %%%%%%%

% The Dirichlet Parameters
alpha = randi(10, 1, k);
% alpha = [0.1, 0.1, 0.2, 0.3, 0.3];

% Topic-topic co-occurrence matrix
[R, alpha] = gen_matrix_R(alpha, 0.3);

% The word-word co-occurrence matrix
Q = A * R * A';
Q_bar = rowStoc(Q);

%% Initialization

% The types of anchors
type_anchors = zeros(n_iter, n_anchor);

% The Q
array_err_Q = zeros(n_iter, 1);
Q_reg = zeros(n_iter, V, V);
array_err_Q_KL = zeros(n_iter, V);

% The A_rec
array_err_A = zeros(n_iter, 1);
A_rec_reg = zeros(n_iter, V, k);

% Best and Worst set of anchors
array_err_A_best = zeros(n_iter, 1);
array_err_A_worst = zeros(n_iter, 1);
A_rec_best_reg = zeros(n_iter, V, k);
A_rec_worst_reg = zeros(n_iter, V, k);

% Without normalizing the docs
if flag_no_norm
    A_rec_best_reg_no_norm = zeros(n_iter, V, k);
end

% Adjusted anchors
if flag_adjust
    type_anchors_adjusted = zeros(n_iter, n_anchor);
    A_rec_reg_adjusted = zeros(n_iter, V, k);
    array_err_A_adjusted = zeros(n_iter, 1);
end

% With true anchors
if flag_trueanchors
    A_rec_reg_trueanchors = zeros(n_iter, V, k);
    array_err_A_trueanchors = zeros(n_iter, 1);
end

%% RUN!!
tic
for i_iter = 1:n_iter
    
    % Display the current iteration
    if mod(i_iter, 10) == 0
        i_iter
        toc
    end
    
    % Topics
    topics = drchrnd(alpha, n_Doc);
    R_empirical = topics' * topics / n_Doc;
    
    % Documents
    x = gen_Docs(topics, A, l_Doc, 'fixed');
    
    %% Analyze the matrix Q
    
    % Word-Word Co-occurrance
    [Q_emp, Q_bar_emp] = gen_matrix_Q(x, 0, 'original');
    Q_reg(i_iter, :, :) = Q_emp;
    
    % Normalized Squared Error
    array_err_Q(i_iter) = norm(Q - Q_emp, 'fro') / norm(Q, 'fro');
    % array_err_Q_KL(i_iter, :) = compute_err(Q_bar, Q_bar_emp, 'row_KL')';
    array_err_Q_l2(i_iter, :) = compute_err(Q_bar, Q_bar_emp, 'row_l2')';    
    %% Recover the topics
    
    % Find the anchors
    candidates = 1:V;
    [anchor_inds, anchors] = find_anchors(Q_emp, candidates, k, 0);
    % [anchor_inds_conv, anchors] = find_anchors_conv(Q_emp, candidates, k, 0);
    
    % Process the anchors
    anchor_inds = sort(anchor_inds);
    for i_k = 1:n_anchor
        type_anchors(i_iter, i_k) = sum(anchor_inds <= k  *i_k);
    end
    for i_k = n_anchor:-1:2
        type_anchors(i_iter, i_k) = type_anchors(i_iter, i_k) - type_anchors(i_iter, i_k-1);
    end
    
    anchor_modk = mod(anchor_inds, k);
    anchor_modk(anchor_modk == 0) = k;
    [~, ind_sort] = sort(anchor_modk);
    anchor_inds = anchor_inds(ind_sort);
    
    % Recover the matrix A
    [A_rec, R_rec] = recoverKL(Q_emp, anchor_inds);
    % [A_rec, R_rec] = recover(Q_emp, anchor_inds)
    A_rec_reg(i_iter, :, :) = A_rec;
    
    array_err_A(i_iter) = norm(A - A_rec, 'fro') / norm(A, 'fro');
    
    %% Recover with better/worst set of anchors
    [A_rec_best, R_rec] = recoverKL(Q_emp, best_anchors);
    array_err_A_best(i_iter) = norm(A - A_rec_best, 'fro') / norm(A, 'fro');
    
    [A_rec_worst, R_rec] = recoverKL(Q_emp, worst_anchors);
    array_err_A_worst(i_iter) = norm(A - A_rec_worst, 'fro') / norm(A, 'fro');
    
    A_rec_best_reg(i_iter, :, :) = A_rec_best;
    A_rec_worst_reg(i_iter, :, :) = A_rec_worst;
    
    %% No Normalization
    if flag_no_norm
        [Q_emp, Q_bar_emp] = gen_matrix_Q(x, 0, 1);
        [A_rec_best, R_rec] = recoverKL(Q_emp, best_anchors);
        array_err_A_best_no_norm(i_iter) = norm(A - A_rec_best, 'fro') / norm(A, 'fro');
        A_rec_best_reg_no_norm(i_iter, :, :) = A_rec_best;
    end
    
    %% Adjusted anchors
    if flag_adjust
        anchor_inds = adjust_anchors(Q, anchor_inds);
        
        % Process the anchors
        anchor_inds = sort(anchor_inds);
        for i_k = 1:n_anchor
            type_anchors_adjusted(i_iter, i_k) = sum(anchor_inds <= k  *i_k);
        end
        for i_k = n_anchor:-1:2
            type_anchors_adjusted(i_iter, i_k) = type_anchors_adjusted(i_iter, i_k) - type_anchors_adjusted(i_iter, i_k-1);
        end
        
        anchor_modk = mod(anchor_inds, k);
        anchor_modk(anchor_modk == 0) = k;
        [~, ind_sort] = sort(anchor_modk);
        anchor_inds = anchor_inds(ind_sort);
        
        [A_rec, R_rec] = recoverKL(Q_emp, anchor_inds);
        A_rec_reg_adjusted(i_iter, :, :) = A_rec;
        array_err_A_adjusted(i_iter) = norm(A - A_rec, 'fro') / norm(A, 'fro');
    end
    
    %% True anchors
    if flag_trueanchors
        true_anchors = Q_bar(best_anchors, :);
        [A_rec_trueanchors, R_rec_trueanchors] = recoverKL(Q_emp, anchor_inds, true_anchors);
        
        A_rec_reg_trueanchors(i_iter, :, :) = A_rec_trueanchors;
        array_err_A_trueanchors(i_iter) = norm(A - A_rec_trueanchors, 'fro') / norm(A, 'fro');
    end
end

%% Mean of results

% The Types
mean_types = mean(type_anchors);

% The Q
mearn_array_err_Q = mean(array_err_Q)
Q_emp_mean = reshape(mean(Q_reg), [V, V]);
err_Q_emp_mean = norm(Q - Q_emp_mean, 'fro') / norm(Q, 'fro')

% The A_rec
mean_err_A = mean(array_err_A)
A_rec_mean = reshape(mean(A_rec_reg), [V, k]);
err_A_rec_mean = norm(A - A_rec_mean, 'fro') / norm(A, 'fro')

% The A_rec_best
mean_err_A_best = mean(array_err_A_best)
A_rec_best_mean = reshape(mean(A_rec_best_reg), [V, k]);
err_A_rec_best_mean = norm(A - A_rec_best_mean, 'fro') / norm(A, 'fro')

% The A_rec_worst
mean_err_A_worst = mean(array_err_A_worst)
A_rec_worst_mean = reshape(mean(A_rec_worst_reg), [V, k]);
err_A_rec_worst_mean = norm(A - A_rec_worst_mean, 'fro') / norm(A, 'fro')

%% Recover using the empirical Q
[A_rec, R_rec] = recoverKL(Q_emp_mean, best_anchors);
err_usingMeanQ = norm(A - A_rec, 'fro') / norm(A, 'fro')
error_between_A = norm(A_rec_best_mean - A_rec, 'fro')./ norm(A, 'fro')

%% Plot the histograms
if flag_histogram
figure;
hist(array_err_A)
hold on;
title('array_{err_A}');

figure;
hist(array_err_A_best)
hold on;
title('array_{err_A,best}');

figure;
hist(array_err_A_worst)
hold on;
title('array_{err_A,worst}');
end

%% Correlations

% err_Q v.s. err_A
figure;
hold on;
title('err_Q v.s. err_A');
scatter(array_err_Q, array_err_A);
scatter(array_err_Q, array_err_A_best, 'k');
scatter(array_err_Q, array_err_A_worst, 'r');
legend('find anchors', 'best anchors', 'worst anchors')
xlabel('err_Q')
ylabel('err_A')

% mean_type v.s. err_A
average_type = mean(type_anchors, 2);
figure;
hold on;
title('Average type v.s. err_A');
scatter(average_type, array_err_A);
xlabel('Average type')
ylabel('err_A');

% # of type 1 v.s. err_A
figure;
hold on;
title('# of type 1 anchors v.s. err_A');
scatter(type_anchors(:, 1), array_err_A);
xlabel('# of type 1 anchors')
ylabel('err_A');

%% Adjusted anchors
if flag_adjust
    
    % The Types
    mean_types = mean(type_anchors)
    mean_types_adjusted = mean(type_anchors_adjusted)
    
    % err_Q v.s. err_A
    figure;
    hold on;
    title('err_Q v.s. err_A');
    scatter(array_err_Q, array_err_A);
    scatter(array_err_Q, array_err_A_adjusted, 'ms');
    scatter(array_err_Q, array_err_A_best, 'k');
    scatter(array_err_Q, array_err_A_worst, 'r');
    legend('find anchors', 'adjusted anchors', 'best anchors', 'worst anchors')
    xlabel('err_Q')
    ylabel('err_A')
    
end
%% Exact anchors
if flag_trueanchors
    mean_err_A_trueanchors = mean(array_err_A_trueanchors)
    
    % err_Q v.s. err_A
    figure;
    hold on;
    title('err_Q v.s. err_A');
    scatter(array_err_Q, array_err_A);
    scatter(array_err_Q, array_err_A_trueanchors, 'ms');
    scatter(array_err_Q, array_err_A_best, 'k');
    scatter(array_err_Q, array_err_A_worst, 'r');
    legend('find anchors', 'exact anchors', 'best anchors', 'worst anchors')
    xlabel('err_Q')
    ylabel('err_A')
end

%% No Normalization
if flag_no_norm
    mean_err_A_best_no_norm = mean(array_err_A_best_no_norm)
    A_rec_best_mean_no_norm = reshape(mean(A_rec_best_reg_no_norm), [V, k]);
    err_A_rec_best_mean_no_norm = norm(A - A_rec_best_mean_no_norm, 'fro') / norm(A, 'fro')
end

%% Row sum of l2 error
mean_l2 = mean(array_err_Q_l2(:, 1:3*k))
figure;
stem(mean_l2)
meanp = A*alpha';
hold on;
stem(meanp(1:3*k)/100, 'k*');
legend('l2 error', 'true probability')

figure;
hold on;
plot(array_err_Q_l2(:, 1:k)', 'rs' );
plot(array_err_Q_l2(:, k+1:2*k)', 'bo' );
plot(array_err_Q_l2(:, 2*k+1:3*k)', 'k*' );

% legend('worst anchors', 'okay anchors', 'best anchors');
xlabel('anchors');
ylabel('l2 error');

%% KL Divergence
mean_KL = mean(array_err_Q_KL(:, 1:3*k))
figure;
stem(mean_KL)

figure;
hold on;
plot(array_err_Q_KL(:, 1:k)', 'rs' );
plot(array_err_Q_KL(:, k+1:2*k)', 'bo' );
plot(array_err_Q_KL(:, 2*k+1:3*k)', 'k*' );

% legend('worst anchors', 'okay anchors', 'best anchors');
xlabel('anchors');
ylabel('KL Divergence');