function A = gen_matrix_A(V, k, p_anchor)
% gen_matrix_A generates topic distribution matrices A with given
% parameters
% V: The number of words
% k: The number of topics
% p_anchor: probability of each anchor word

% The Autism NMF Project
% Hongyao Ma
% Created:  01/11/2015
% Updated:  01/11/2015

% Summation of anchors
sum_p_anchor = sum(p_anchor);
if sum_p_anchor > 1
    error('Summation of anchors should be no greater than 1!');
end

% Number of anchors
n_anchor = length(p_anchor)*k;
if n_anchor > V
    error('Total number of anchors should be smaller than V!');
end

% The anchors
A1 = [];
for i = 1:length(p_anchor) 
    A1 = [A1; eye(k) * p_anchor(i)];
end

% The rest of the topic distributions
A2 = rand(V - n_anchor, k);
%A2 = randi(2, V-k, k)-1
sumA2 = sum(A2);
A2 = (1 - sum_p_anchor)* A2./ repmat(sumA2, V - n_anchor, 1);

% Concatenate the matrices
A = [A1; A2];

% TODO:
% 1. Different probabilities for different anchors
% 2. Get more than 1 anchor for each topic