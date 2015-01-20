function [A, V, k] = gen_trickyA(A_type)
% gen_trickyA generates toy topic distributions that are difficult to
% distinguish for topic modeling algorithms

% The Autism NMF Project
% Hongyao Ma
% 01/19/2015

%% The tricky topic distributions
% Anchors: 1, 2, 3
A_1 = repmat(eye(3), 1, 3);

% Anchors: 1, 4, 7
A_2 = [1 1 1 0 0 0 0 1 1; ...
       0 1 1 1 1 1 0 0 0; ...
       0 0 0 0 1 1 1 1 1];

% Anchors: 1, 4, 7
A_3 = [2 1 1 0 1 1 0 0 0; ...
       0 0 0 2 1 1 0 1 1; ...
       0 1 1 0 0 0 2 1 1 ];

% No clear anchors
A_4 = [1 1 1 1 1 1 0 0 0; ...
       1 1 1 0 0 0 1 1 1; ...
       0 0 0 1 1 1 1 1 1];

A_5 = [eye(3), ones(3, 6)];

A_6 = [eye(3)*3, eye(3)*2, eye(3)];
A_7 = [eye(3)*9, eye(3)*3, eye(3)];

%% Choose the one specified
switch A_type
    case 1
        A = A_1';
    case 2
        A = A_2';
    case 3
        A = A_3';
    case 4
        A = A_4';
    case 5
        A = A_5';
    case 6
        A = A_6';
    case 7
        A = A_7';        
    otherwise
        error('Type of A not found!');
end

%% Parameters and Normalization
% Sizes
V = size(A, 1);
k = size(A, 2);

% Normalize
A = A ./ repmat(sum(A), V, 1);