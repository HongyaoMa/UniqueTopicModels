function err_out = compute_err(x_true, x_emp, err_type)
% compute_err compute the error of x_emp compared with x_true, of
% corresponding types specified by err_type

% The Autism NMF Project
% Hongyao Ma
% Created:  01/16/2015
% Modified: 01/16/2015

% Default error type
if nargin == 2
    display('Error Type Not Specified! Using L2 error by default!');
    err_type = 'l2';
end

switch err_type
    % The L1 error
    case 'l1'       
        err_out = abs(x_true-x_emp);
    
    % The l2 error
    case 'l2'
        err_out = (x_true - x_emp).^2;
        
    % The rwo sum of l2 error
    case 'row_l2'
        err_out = sum((x_true - x_emp).^2, 2);        
    
    % The KL divergence for each row 
    case 'row_KL'
        rowsum = [sum(x_true, 2); sum(x_emp, 2)];
        if ~all(abs(rowsum - 1) < 1e-2)
            error('Rows do not sum to 1');
        end
        if any(any ((x_true == 0) .* (x_emp ~= 0)))
            error('KL Divergence Not Defined!')
        end
        err_out = sum(x_true .* log(x_true ./ x_emp), 2);
    otherwise
        error('Unspecified error type!!');
end

end

%% TODO: 
% 1. Better error message
% 2. Add the option to normalize w.r.t. the true distribution

