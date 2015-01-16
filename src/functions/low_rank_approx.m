function [array_err] = low_rank_approx(x, ranks, plot_error)
% err_low_r_approx computes and plots the error of the low-rank
% approximations

% The Autism NMF Project
% Hongyao Ma
% Created:  01/11/2015
% Modified: 01/11/2015

% The maximum rank to consider
if nargin == 1
    ranks = 1 : min(size(x));
end

% Check that the max_r is within the size of the input
if max(ranks) > min(size(x))
    display('Warning: Specified maximum rank is too large!');
    display('Truncated at the size of the input matrix.');
    ranks = ranks(ranks <= min(size(x)));
end

cutfirst = 0;
if ranks(1) ~= 1
    ranks = [1, ranks];
    cutfirst = 1;
end

% SVD 
[U, S, V] = svd(x);
array_err = zeros(size(ranks));

% Rank-1 approximation
x_approx = U(:, 1) * S(1,1) * V(:,1)';
array_err(1) = norm(x - x_approx, 'fro');

% Compute the error of low-rank approximations
for i_rank = 1:length(ranks)-1
    x_approx = x_approx + U(:, ranks(i_rank)+1 : ranks(i_rank+1)) ...
        * S(ranks(i_rank) + 1 :  ranks(i_rank + 1), ranks(i_rank) + 1: ranks(i_rank + 1))...
        * V(:, ranks(i_rank) + 1:  ranks(i_rank + 1))';
    array_err(i_rank + 1) = norm(x - x_approx, 'fro');
end

% Normalize by the frobenius norm of the input
array_err = array_err / norm(x, 'fro');

if cutfirst
    array_err = array_err(2:end);
    ranks = ranks(2:end);
end

% Plot the results
if plot_error
    figure;
    plot(ranks, array_err, 'o--');
    xlabel('rank');
    ylabel('err of low-rank approx');
end


