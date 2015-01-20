function alpha = gen_alpha(alpha0, k, type)
% gen_alpha generates dirichlet parameters of given type

% The Autism NMF Project
% Hongyao Ma
% 01/19/2015

switch type
    case 'uniform'
        alpha = ones(1, k);
    case 'random'
        alpha = randi(10, 1, k);
    otherwise
        error('Undefined type of Dirichlet parameters!');
end

% Normalization
alpha = alpha / sum(alpha) * alpha0;

end