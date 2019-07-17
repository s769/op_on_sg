format long


% 1. Coefficients of polynomials with respect to lambda

% Specify the dimension of the polynomial we're looking at
% For some reason this is actually coefficients of polynomial 
% of dimension (dim-1)
dim = 6;

% Generate input values for lambda
% lambda is chosen from 10^(-12) to 10^(12) for all integer 
% powers of 10
lambda = zeros(1, 25);
for i=1:1:25
    lambda(i) = 10^(i-13);
end

% Initialize space to store the coefficients of the matrix
% Matrix has dimension dim*25, i-th row represent coefficient of P_{ij}
store_mat = zeros(dim, 25);


% Loop to generate values for different r
for i = 1:1:25
    % Due to current instability in code, we need to generate more
    % rows than needed
    sobolov_vals = sobolov_coeff(dim+3, lambda(i), 3, 0);
    % Take the right terms in the matrix
    store_mat(:, i) = sobolov_vals(dim, 1:dim)' ;
end


% Plot the values [log(x)-log(y) plot]
% For plotting a single row
scatter(log10(lambda), log10(abs(store_mat(1,:))));

% For plotting multiple rows
for i=1:1:(dim-1)
    scatter(log10(lambda), log10(abs(store_mat(i,:))))
    hold on
    xlabel('log10(lambda)')
    ylabel('log10(P)')
end

% Finding slope of a row
slope = (log10(abs(store_mat(1, end))) - log10(abs(store_mat(1, 10)))) / (25 - 10)


%{
% 2. Coefficients of polynomials with respect to layer
dim = 6;
lambda = 1;

sobolov_vals = sobolov_coeff(dim+3, lambda, 3, 1);
sobolov_vals = sobolov_vals(1:dim, 1:dim)
%}
