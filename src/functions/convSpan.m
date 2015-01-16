function [ind_covered, num_covered] = convSpan(data, basis)
% convSpan computes the linear span coefficient of basis for each vector in
% data, then examine how many of them have all non-negative coefficients
% When the basis are redundant, the method fails

% The NMF Autism Project
% Hongyao Ma
% 12/11/2014


size(basis)
size(data)

coeff = basis \ data;
size(coeff)
ind_covered = all(coeff > -1e-8); 

num_covered = sum(ind_covered);

coeff_nonneg = coeff;
coeff_nonneg(coeff < 0) = 0;


