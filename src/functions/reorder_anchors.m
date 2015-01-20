function anchor_inds_out = reorder_anchors(anchors_in)
% reorder_anchors reorders anchors w.r.t. the residue after mod k where k
% is the number of anchors
% this helps reordering the columns of the topic distribution matrix so
% that the result can be compared easily with the ground truth

% The Autism NMF Project
% Hongyao Ma
% 01/20/2015

% Number of anchors
k = length(anchors_in);

% Mod k
anchor_modk = mod(anchors_in, k);

% Change the "0" to k
anchor_modk(anchor_modk == 0) = k;

% Reorder the anchors 
[~, ind_sort] = sort(anchor_modk);
anchor_inds_out = anchors_in(ind_sort);

end