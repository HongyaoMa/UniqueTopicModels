function scatter3D(x)
% scatter3D plots a N by 3 matrix as a 3D scatter plot

% The Autism NMF Project
% Hongyao Ma
% 10/13/2014

if size(x,1) == 3
    x = transpose(x);
end

if size(x, 2) ~= 3
    error('Input should be N by 3 matrices')
end

figure;
scatter3(x(:,1), x(:,2), x(:,3), '.');