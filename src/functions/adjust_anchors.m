function anchor_inds_out = adjust_anchors(Q, anchor_inds_in, range)
% adjust_anchors finds alternative sets of anchors that are close to the
% original ones but has higher probabilities

% The Autism NMF Project
% Hongyao Ma
% Created:  01/14/2015
% Modified: 01/14/2015

if nargin == 2
    range = 0.01;
    display('Range not specified. Set to default = 0.01');
end

% Normalize to get the matrix Q_bar
Q_rowsum = sum(Q, 2);
Q_bar = Q ./ repmat(Q_rowsum, 1, size(Q, 2));

% Number of anchors
k = length(anchor_inds_in);

% Initialize the output indices
anchor_inds_out = anchor_inds_in;

for i_anchor = 1:k
    
    % The current anchor
    curr_anchor = Q_bar(anchor_inds_in(i_anchor), :);
    
    % Find the alternatives within the range
    Q_diff = Q_bar - repmat(curr_anchor, length(Q), 1);
    dist = sum(Q_diff.^2, 2) ./ sum(curr_anchor.^2);
    ind_close_words = find(dist < range);
    
    % Probability of the alternatives
    p_close_words = Q_rowsum(ind_close_words);

    % Pick the one with highest probability as the new anchor
    [~, ind_close_words_maxp] = max(p_close_words);    
    anchor_inds_out(i_anchor) = ind_close_words(ind_close_words_maxp);
    
end

end