function [anchor_indices, anchors] = find_anchors(Q, candidates, n_anchor, plot_result)
% find_anchors implements the anchor word finding algorithm in the topic
% modeling paper, which was previously implemented in Python

% The Autism NMF Project
% Hongyao Ma
% Created:  12/11/2014
% Modified: 12/15/2014

% The matrix Q_bar
Q_rowsum = sum(Q, 2);
Q_bar = Q ./ repmat(Q_rowsum, 1, size(Q, 2));

% Take only the candidates
Q_bar_can = Q_bar(candidates, :);

% The original candidates
Q_bar_can_orig = Q_bar_can;

d = size(Q_bar, 2);
n_words = length(candidates);

% Initialization
anchors = zeros(n_anchor, d);
anchor_indices = zeros(n_anchor, 1);
basis = zeros(n_anchor-1, d);
mean_dist = zeros(n_anchor, 1);
max_dist = zeros(n_anchor, 1);

%% 1st Anchor
% Compute the distances from the origin
dist = sum(Q_bar_can.^2, 2);
[maxDist, ind_max] = max(dist);

mean_dist(1) = mean(dist);
max_dist(1) = maxDist;

% Take the first basis vector as the one farthest away from the origin
anchors(1, :) = Q_bar_can(ind_max, :);
anchor_indices(1) = ind_max;

% Make anchor 1 the origin of the coordinate system
Q_bar_can = Q_bar_can - repmat(anchors(1,:), n_words, 1);

% Plot the current Q matrix and anchor
if plot_result
    figure;
    subplot(1,2,2);        hold on;
    scatter3(Q_bar_can(:, 1), Q_bar_can(:, 2), Q_bar_can(:, 3), 'b*');
    scatter3(Q_bar_can(anchor_indices(1), 1), Q_bar_can(anchor_indices(1), 2), Q_bar_can(anchor_indices(1), 3), 'ro');
    subplot(1,2,1);        hold on;
    scatter3(Q_bar(:, 1), Q_bar(:, 2), Q_bar(:, 3), 'b*');
    scatter3(Q_bar_can_orig(anchor_indices(1), 1), Q_bar_can_orig(anchor_indices(1), 2), Q_bar_can_orig(anchor_indices(1), 3), 'ro');
end

%% 2nd Anchor
% Take the 2nd basis vector as the one farthest from the first one
dist = sum(Q_bar_can.^2, 2);
[maxDist, ind_max] = max(dist);
anchors(2, :) = Q_bar_can_orig(ind_max, :);
anchor_indices(2) = ind_max;

mean_dist(2) = mean(dist);
max_dist(2) = maxDist;

% First basis vecror
basis(1,:) = Q_bar_can(anchor_indices(2), :) / norm(Q_bar_can(anchor_indices(2), :));

% Find the rest of the anchors
for i = 2:n_anchor-1
    coeffs = Q_bar_can * basis(i-1, :)';
    Q_bar_can = Q_bar_can - diag(coeffs) * repmat(basis(i-1, :), n_words, 1);
    
    % Plot the current result
    if plot_result
        figure;
        subplot(1,2,1);        hold on;
        scatter3(Q_bar(:, 1), Q_bar(:, 2), Q_bar(:, 3), 'b*');
        scatter3(Q_bar(anchor_indices(1:i), 1), Q_bar(anchor_indices(1:i), 2), Q_bar(anchor_indices(1:i), 3), 'ro');
        subplot(1,2,2);        hold on;
        scatter3(Q_bar_can(:, 1), Q_bar_can(:, 2), Q_bar_can(:, 3), 'b*');
        scatter3(Q_bar_can(anchor_indices(1:i), 1), Q_bar_can(anchor_indices(1:i), 2), Q_bar_can(anchor_indices(1:i), 3), 'ro');
    end
    
    if norm(Q_bar_can, 'fro') / norm (Q_bar_can_orig, 'fro') < 1e-7
        display('Rank of input is smaller than number of anchors asked!');
        break
    end
    
    dist = sum(Q_bar_can.^2, 2);
    [maxDist, ind_max] = max(dist);
    anchors(i+1, :) = Q_bar_can_orig(ind_max, :);
    anchor_indices(i+1) = ind_max;
    basis(i,:) = Q_bar_can(anchor_indices(i+1), :) / norm(Q_bar_can(anchor_indices(i+1), :));
    
    mean_dist(i+1) = mean(dist);
    max_dist(i+1) = maxDist;
end

if plot_result
    figure;
    stem(mean_dist);
    hold on;
    stem(max_dist, 'k*');
    legend('mean distance', 'max distance');
end

% TODO: Adjust the indices when the candidates are not all

