function [sob_coeff_mat] = sobolov_coeff(r, lambda, j, norm_flag)
%SOBOLOV_COEFF 
%   Goal: Gives a representation of Sobolev polynomials in terms
%           of the monomial basis
%   Inputs: r is the number of orders we want to calculate the coefficients
%           lambda is the coefficient used in the Sobolev inner product
%           j\in {1,2,3} is the type of SOP that we're working with
%           norm_flag = 1 if we want to normalize, 0 if not
%   Outputs: sob_coeff_mat is a (r+1)*(r+1) matrix representing the
%               coefficients of the Sobolev orthogonal polynomials               


% First build up the inner product matrix
orig_inner_val = inner_prod_orig(r+1, j, j, 1, 1, 1);
sob_inner_val = inner_prod_sob(r+1, j, j, 1, 1, 1, 1);


% Allocate space for the final matrix
sob_coeff_mat = zeros(r+1, r+1);

% Main loop for generating the polynomials
for l1=1:1:(r+1)
    for l2=1:1:(r+1)
        % Make all 1's on the diagonal [since we're dealing with monic SOP]
        if (l1 == l2)
            sob_coeff_mat(l1, l2) = 1;
        end
        
        % Main recursion part
        if (l1 > l2)
            % First recurse on the outer loop sum
            for l3=(l2-1):1:(l1-2)
                % First find the numerator
                temp_num = 0;
                
                % Part I
                for t=0:1:l3
                    temp_num = temp_num + sob_coeff_mat(l3+1, t+1) ...
                        * orig_inner_val(l1, t+1);
                end
                
                % Part II
                for t=1:1:l3
                    temp_num = temp_num + sob_coeff_mat(l3+1, t+1) ...
                        * lambda * orig_inner_val(l1-1, t);
                end
                
                % Then find the denominator
                temp_denom = 0;
                % Part I
                for t1=0:1:l3
                    for t2=0:1:l3
                        temp_denom = temp_denom + sob_coeff_mat(l3+1, t1+1) ...
                            * sob_coeff_mat(l3+1, t2+1) * ...
                            orig_inner_val(t1+1, t2+1);
                    end
                end
                
                % Part II
                for t1=1:1:l3
                    for t2=1:1:l3
                        temp_denom = temp_denom + sob_coeff_mat(l3+1, t1+1) ...
                            * sob_coeff_mat(l3+1, t2+1) * ...
                            orig_inner_val(t1, t2) * lambda;
                    end
                end
                
                % Adding everything together
                sob_coeff_mat(l1,l2) = sob_coeff_mat(l1,l2) + temp_num ...
                    * sob_coeff_mat(l3+1, l2) / temp_denom;
            end
            sob_coeff_mat(l1,l2) = -sob_coeff_mat(l1, l2);
        end
    end
end

% Normalizing the vector if it is requested
if (norm_flag == 1)
    for l1=1:1:r
        % Finding the inner product
        temp_inner = 0;
        % Part I: Ordinary inner product
        for l2=1:1:l1
            for l3=1:1:l1
                temp_inner = temp_inner + sob_coeff_mat(l1, l2) * ...
                    sob_coeff_mat(l1, l3) * orig_inner_val(l2, l3);
            end
        end
        
        % Part II: Sobolev addition
        for l2=2:1:l1
            for l3=2:1:l1
                temp_inner = temp_inner + lambda * sob_coeff_mat(l1, l2) * ...
                    sob_coeff_mat(l1, l3) * orig_inner_val(l2-1, l3-1);
            end
        end
        sob_coeff_mat(l1,:) = sob_coeff_mat(l1,:) / sqrt(temp_inner);
    end
end
 

end

