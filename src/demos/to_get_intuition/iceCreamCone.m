% Demo Ice Cream Cone

% Non-negative Matrix Factorization
% Hongyao Ma
% 10/01/2014

clear;
clc;
close all;

% Generate Random Points within [0, 1]^p
p = 3;
N = 2000;
x = rand(N, p);
% x = [1 1 1]
x(:,3) = 1;

% Get the ice-cream cone
ind_iceCreamCone = sum(x, 2) >= (p-1)^0.5 * sum(x.^2, 2).^0.5;
sum(ind_iceCreamCone)
iceCreamCone = x(ind_iceCreamCone, :);

% Plot the cone
figure;
if p == 2
    scatter(iceCreamCone(:,1), iceCreamCone(:,2), '.');
elseif p == 3
    scatter3(iceCreamCone(:,1), iceCreamCone(:, 2), iceCreamCone(:,3), '.');
end

