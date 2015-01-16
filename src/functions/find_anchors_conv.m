function [anchor_indices, anchors, rankplus, more_anchors, more_inds] = find_anchors_conv(Q, candidates, n_anchor, plot_result)
% find_anchors_conv finds the "anchor words" that are approximately the
% vertices of the convex hulls spanned by the input

% The Autism NMF Project
% Hongyao Ma
% Created: 12/15/2014
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
% Q_bar_can = Q_bar_can - repmat(anchors(1,:), n_words, 1)

% Plot the current Q matrix and anchor
if plot_result
    figure;
    hold on;
    scatter3(Q_bar(:, 1), Q_bar(:, 2), Q_bar(:, 3), 'b*');
    scatter3(Q_bar_can_orig(anchor_indices(1), 1), Q_bar_can_orig(anchor_indices(1), 2), Q_bar_can_orig(anchor_indices(1), 3), 'ro');
end

%% Rest of the anchors
for i = 2:n_anchor
    
    % Compute the distances
    [ind_in, res_dis, coeffs] = inConvHull(anchors(1:i-1, :)', Q_bar_can');
    [maxDist, ind_max] = max(res_dis);
    
    if maxDist < 1e-7
        display('Rank of input is smaller than number of anchors asked!');
        break
    end
    
    % Actual non-negative rank
    rankplus = i;
    
    % Take the one furthest away
    anchors(i, :) = Q_bar_can(ind_max, :);
    anchor_indices(i) = ind_max;
    
    % Plot the current result, if required
    if plot_result
        figure;
        hold on;
        scatter3(Q_bar(:, 1), Q_bar(:, 2), Q_bar(:, 3), 'b*');
        scatter3(Q_bar_can_orig(anchor_indices(1:i), 1), Q_bar_can_orig(anchor_indices(1:i), 2), Q_bar_can_orig(anchor_indices(1:i), 3), 'ro');
    end
    
end

rankplus = n_anchor;

% Comute the upper bound of the non-negative rank of the matrix
% Looking for rankplus anchors that convexly span the data
% Not how rank+ is defined

% Check how far are datapoints from the convex hull
[~, res_dis] = inConvHull(anchors', Q_bar_can');
[maxDist, ind_max] = max(res_dis);

more_anchors = anchors;
more_inds = anchor_indices;

% Keep looking until exhaust the data
while maxDist > 1e-7
    
    % Increase rankplus and record the newly found anchor
    rankplus = rankplus + 1;
    more_anchors = [more_anchors ; Q_bar_can(ind_max, :)]; 
    more_inds = [anchor_indices; ind_max];

    % Update the residual distances
    [~, res_dis] = inConvHull(more_anchors', Q_bar_can');
    [maxDist, ind_max] = max(res_dis);
end

% TODO: adjust the indices
    