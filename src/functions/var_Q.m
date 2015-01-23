function var_analytical = var_Q(p, l_Doc)

% Total number of words
V = length(p);

pi2pj2 = p.^2 * (p.^2)';
pi2pj = p.^2 * p';
pipj2 = p * p.^2';
pipj = p*p';
pi4 = p.^4;
pi3 = p.^3;
pi2 = p.^2;

% The variances
var_analytical = (-4 * pi2pj2 + pi2pj + pipj2 ) * l_Doc^3 + ...
    (10 * pi2pj2 - 3* pi2pj - 3* pipj2 + pipj) * l_Doc^2 + ...
    (-6 * pi2pj2 + 2* pi2pj + 2* pipj2 - pipj) * l_Doc;

% Diagonal elements
var_analytical_diag = l_Doc * (-2 *pi2 +8 * pi3 -6 * pi4) +...
    l_Doc^2 * (2 * pi2 - 12* pi3 + 10 * pi4) + ...
    l_Doc^3 * (4 * pi3 - 4 * pi4);
var_analytical(1:V+1:end) = var_analytical_diag;

% Normalized each document by l_doc(l_doc-1)
var_analytical = var_analytical / l_Doc^2 / (l_Doc - 1)^2;

%% 
Epi2 = 


Epi2pj2 = 