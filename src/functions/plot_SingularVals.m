function plot_SingularVals(x, max_r)
% plot_SingularVals plots the singular values of the matrix x

% The Autism NMF Project
% Hongyao Ma
% 01/11/2015

% The maximum rank to consider
if nargin == 1
    max_r = min(size(x));
end

% Check that the max_r is within the size of the input
if max_r > min(size(x))
    error('Specified maximum rank is too large!');
end

% Compute the singuolar values
s = svds(x, max_r);

% Plot in log scale if the difference is really big
figure;
if max(s) / min(s) > 1000
    semilogy(1:length(s), s, 'o');
else
    stem(s);
end







