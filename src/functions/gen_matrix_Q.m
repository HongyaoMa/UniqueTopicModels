function [Q, Q_bar] = gen_matrix_Q(data, show_results, type_normalization)
% gen_matrix_Q generates the co-occurrance matrix Q from input data

% The Autism NMF Project
% Hongyao Ma
% Created:  11/25/2014
% Updated:  01/14/2015

% Set the default to "normalize the docs"
if nargin == 2
    type_normalization = 'original';
end

% Size
[m, n] = size(data);

% Words per document
wpd = sum(data, 2);

switch type_normalization
    case 'original'
        norm_mat = diag((1./wpd./(wpd-1)).^0.5);
    case 'none'
        norm_mat = eye(m);
    case 'n32'
        norm_mat = diag(wpd.^(-3/4));
    case 'n'
        norm_mat = diag(wpd.^(-0.5));        
    otherwise
        error('Normalization Type Unkonwn!');
end

data_2_norm = norm_mat * data;

%% Q - the co-occurrance matrix
Q = data_2_norm' * data_2_norm - diag(sum(norm_mat * data_2_norm));
Q = Q/(sum(sum(Q)));

% Visualize the Q
if show_results
    figure;
    imshow(Q * m); % The numbers are too small
    title('Q');
end

%% Q_bar

% Normalize the rows
Q_row_sum = 1 ./ sum(Q, 2);
Q_bar = Q .* repmat(Q_row_sum, 1, size(Q,2));

% Visualize the Q
if show_results
    figure;
    imshow(Q_bar * n/2);
    title('Q_bar');
    
end
