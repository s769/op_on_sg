function [sob_inner_val] = inner_prod_sob(r, m, n, s, t, lambda, mat_flag)
%INNER_PROD_VAL 
%   Goal: Generate the inner product matrix of Sobolov norm
%       A_{mn}(s+1,t+1) = <P_sm, P_tn>+lambda<P_(s-1)m, P_(t-1)n>
%   Inputs: r represent the number of layers we want to calculate
%           1<=m,n<=3 represent the second index we're working on
%           0<=s,t<=r represent the first index we're working on
%           lambda represent the coefficient in the Sobolov norm
%           mat_flag is 1 when we want to find the entire matrix, 0 if only
%               1 value is needed
%   Outputs: inner_val is an (r+1)*(r+1) matrix representing all the inner
%               products, or a single value if only one is needed

% Fetch the original inner products
orig_inner_val = inner_prod_orig(r, m, n, s, t, 1);

% Initialize the matrix
sob_inner_val = orig_inner_val;

% Main loop of iterating over all values
for j = 2:1:(r+1)
    for k = 2:1:(r+1)
        sob_inner_val(j, k) = sob_inner_val(j, k) + lambda * ...
            sob_inner_val(j-1, k-1);
    end
end

if (~mat_flag)
    sob_inner_val = sob_inner_val(s+1, t+1);
end

end

