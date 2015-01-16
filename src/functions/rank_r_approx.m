function err = rank_r_approx(x, r)
% rank_r_approx computes and error of the rank-r approximation to the
% imput matrix x

% The Autism NMF Project
% Hongyao Ma
% Created:  01/11/2015
% Modified: 01/11/2015

if r > min(size(x))
    error('Rank is too large!');
end

[U, S, V] = svd(x);
x_r = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';

err = norm(x-x_r, 'fro') / norm(x, 'fro');