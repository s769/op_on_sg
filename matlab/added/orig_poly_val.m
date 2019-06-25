function [alpha, beta, gamma, eta] = orig_poly_val(r, type)
%ORIG_POLY_VAL Summary of this function goes here
%   Goal: Calculate the values of P_{j1}(q_1), P_{j2}(q_1), P_{j3}(q_1), 
%       \partial_n P_{j1}(q_1) as to generate all values
%   Inputs: r is the number of orders we want to calculate the coefficients
%           type is true(1) when entire list is needed, or 0 when only one
%                value is needed
%   Outputs: alpha is a 1*(r+1) vector denoting P_{j1}(q1) for 0<=j<=r
%            beta is a 1*(r+1) vector denoting P_{j2}(q1) for 0<=j<=r
%            gamma is a 1*(r+1) vector denoting P_{j3}(q1) for 0<=j<=r
%            eta is a 1*(r+1) vector denoting \partial_n P_{j1}(q1) for 0<=j<=r
%            [if for number version, only one number for each output]
%            [For indices, alpha(j)=alpha[j-1]]


% Initialize space to store the variables
% Since gamma depend on the previous value of alpha, to get gamma
%   we need to have an additional alpha calculated.
alpha = zeros(1, r+2);
beta = zeros(1, r+1);
gamma = zeros(1, r+1);
eta = zeros(1, r+1);

% Initialize values of the 
alpha(1) = 1;
alpha(2) = 1/6;
beta(1) = -1/2;
gamma(1) = 3 * alpha(2);
eta(1) = 0;


for j=2:1:(r+1)
    % Update alpha(j+1)
    for l=1:1:(j-1)
        alpha(j+1) = alpha(j+1) + alpha(j-l+1) * alpha(l+1);  
    end 
    alpha(j+1) = (alpha(j+1) * 4) / (5^j-5);
    
    % Update beta(j)
    for l=0:1:(j-2)
        beta(j)=beta(j)+alpha(j-l)*beta(l+1)*(3*5^(j-1-l)-5^(l+1)+6);
    end
    beta(j) = (beta(j) * 2) / (15 * (5^(j-1) - 1));
    
    % Update gamma(j)
    gamma(j) = 3 * alpha(j+1);
    
    % Update eta(j)
    for l=0:1:(j-2)
        eta(j) = eta(j) + eta(l+1) * beta(j-l);
    end
    eta(j) = 2 * eta(j);
    eta(j) = eta(j) + (5^(j-1)+1)*alpha(j) / 2;
end

% Reduce alpha to proper dimension
alpha = alpha(1:r+1);

% Case that we only need values from one iteration
if type == false
    alpha = alpha(r+1);
    beta = beta(r+1);
    gamma = gamma(r+1);
    eta = eta(r+1);
end
end

