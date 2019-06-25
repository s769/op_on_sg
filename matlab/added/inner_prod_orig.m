function [inner_val] = inner_prod_orig(r, m, n, s, t, mat_flag)
%INNER_PROD_VAL 
%   Goal: Generate the inner product matrix of ordinary norm
%       A_{mn}(s+1,t+1) = <P_sm, P_tn>
%   Inputs: r represent the number of layers we want to calculate
%           1<=m,n<=3 represent the second index we're working on
%           0<=s,t<=r represent the first index we're working on
%           mat_flag is 1 when we want to find the entire matrix, 0 if only
%               1 value is needed
%   Outputs: inner_val is an (r+1)*(r+1) matrix representing all the inner
%               products, or a single value if only one is needed

% Store the original polynomial value arrays
[alpha, beta, ~, eta] = orig_poly_val(2*r+3, 1);

% Create alpha_prime
alpha_prime = alpha;
alpha_prime(1) = 1/2;

% Initialize the matrix
inner_val = zeros(r+1, r+1);

% Main loop of iterating over all values
for j = 1:1:(r+1)
    for k = 1:1:(r+1)
        % Find the lower order index
        lower_ind = j - min(j, k);
        
        % Case 1: <P_j1, P_k1>
        if (m == 1 && n == 1)
            for l = lower_ind:1:(j-1)
                inner_val(j,k) = inner_val(j, k)+alpha(j-l)*eta(k+l+1) ...
                    -alpha(k+l+1)*eta(j-l);
            end
            inner_val(j, k) = inner_val(j, k) * 2;
        end   
        
        % Case 2: <P_j2, P_k2>
        if (m == 2 && n == 2)
            for l = lower_ind:1:(j-1)
                inner_val(j,k) = inner_val(j, k)-beta(j-l)*alpha(k+l+1) ...
                    +beta(k+l+1)*alpha_prime(j-l);
            end
            inner_val(j, k) = inner_val(j, k) * 2;
        end
        
        % Case 3: <P_j3, P_k3>
        if (m == 3 && n == 3)
            for l = lower_ind:1:(j-1)
                inner_val(j,k) = inner_val(j, k)+alpha(j-l+1)*eta(k+l+2) ...
                    -alpha(k+l+2)*eta(j-l+1);
            end
            inner_val(j, k) = inner_val(j, k) * 18;
        end
        
        % Case 4: <P_j1, P_k2>
        if (m == 1 && n == 2)
            for l = 0:1:(j-1)
                inner_val(j,k) = inner_val(j, k)-alpha(j-l)*alpha(k+l+1) ...
                    -beta(k+l+1)*eta(j-l);
            end
            inner_val(j, k) = inner_val(j, k) * 2;
        end
        
        % Case 5: <P_j2, P_k1>
        if (m == 2 && n == 1)
            for l = 0:1:(j-1)
                inner_val(k,j) = inner_val(k,j)-alpha(j-l)*alpha(k+l+1) ...
                    +beta(k+l+2)*eta(j-l+1);
            end
            inner_val(k, j) = inner_val(k, j) * 2;
        end
    end
end

if (~mat_flag)
    inner_val = inner_val(s+1, t+1);
end

end

