function [outputArg1,outputArg2] = gram_schmidt(r, lambda, ...
    j, norm_type, norm_flag)
%GRAM_SCHMIDT Implements the Gram-Schmidt algorithm
%   Detailed explanation goes here
%   Inputs: r is the dimension of the matrix we're working with
%           lambda is the coefficient used in the Sobolov inner product
%           j\in {1,2,3} is the type of SOP that we're working with
%           norm_type is 1 for ordinary norm, and 0 for Sobolov
%           norm_flag is 1 for normalization of rows, 0 for not

% Start with the identity matrix
coeff_mat = eye(r+1);

% Find the inner product matrix
if (norm_type == 1)
    inner_mat = inner_prod_orig(r+1, j, j, 1, 1, 1);
else
    inner_mat = inner_prod_sob(r+1, j, j, 1, 1, lambda, 1);
end

% Perform the Gram-Schmidt update
for l1=2:1:(r+1)
    for l2=1:1:(r+1)
        
    end
end

