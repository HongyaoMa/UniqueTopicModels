function x = linear_comb(basis, coeff)
% Compute the linear combinations of columns of "basis" with coefficients
% "coeff"

% The Autism NMF Project
% Hongyao Ma
% Created:  10/13/2014
% Modified: 12/11/2014

% d = size(basis, 1);
r = size(basis, 2);

rr = size(coeff, 1);
% N = size(coeff, 2);

% Check that the sizes match
if r ~= rr
    error('Dimensionality does not match!');
end

% Combinations of columns of basis
x = basis * coeff;

end
