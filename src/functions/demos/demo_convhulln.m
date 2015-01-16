% Demo of convhulln
% Plot for the 2D and 3D cases
% Tried out higher-dimensional cases but cannot plot

% The Autism-NMF Project
% Hongyao Ma
% 10/13/2014

clear;
clc;
close all;

% Add all subfolders of the source directory
addpath(genpath('../src'));


% Generate Random Points within [0, 1]^p
p = 3;
N = 2000;
x = rand(p, N);
% scatter3D(x);

%
X = x';
X = X(1:100, :);
K = convhulln(X, {'Qt'});

%K = convhull(X(:,1), X(:,2));

if p == 2
    plot(X(K, 1),X(K, 2),'k-o',X(:,1),X(:,2),'b.') 
    
elseif p == 3
    figure;
    trisurf(K, X(:,1),X(:,2),X(:,3))

    hold on;
    scatter3(X(:,1), X(:,2), X(:,3), 'ro');    
end
