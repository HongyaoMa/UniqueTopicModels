function [ind_in, res_norm, coeffs] = inConvHull(basis, data)
% inConvHull conputes the distance between the vectors in data and the 
% convex hull spanned by the vectors in basis using linear lease-squares
% method

% The Autism NMF Project
% Hongyao Ma
% 12/11/2014

% Check the input sizes
if size(basis, 1) ~= size(data, 1)
    error('Input dimensions do not match!');
end

% dimension of the data vectors
d = size(basis, 1);

% # of basis vectors
d_basis = size(basis, 2);

% # of data vectors
N = size(data, 2);

% Formulation of linear least-squares problem
C = basis;
A = [];
b = [];

Aeq = ones(1, d_basis);
beq = 1;
lb = zeros(d_basis, 1);
ub = ones(d_basis, 1);

% Initialize
res_norm = zeros(N, 1);
coeffs = zeros(d_basis, N);

curr_options = optimoptions('lsqlin','Display','off', 'Algorithm', 'active-set');
x0 = [];

% Analyze each data vector in data
for i = 1:N
    d = data(:, i);
    [x, resnorm] = lsqlin(C, d, A, b, Aeq, beq, lb, ub, x0, curr_options);
    %[x, resnorm] = lsqlin(C, d, A, b, Aeq, beq, lb, ub);
    res_norm(i) = sqrt(resnorm);
    coeffs(:, i) = x;

end

% Indices for insiders
ind_in = res_norm <= 1e-7;






