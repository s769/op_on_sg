function y = qcontract(x,q)
% This function takes in a column of coordinate pairs and 
% computes their image under the similarity which contracts
% distance to the point q by a factor of 1/2.
%
% support functions: fi.m
%
% INPUT:
% x - 3^(m+1) by 2 matrix holding x and y coordinates
% n - the number of points
%
% OUTPUT: 
% y - coordinates of the n images
%
% from http://www.math.cornell.edu/~mhall/SGprograms/qcontract.m

n = max(size(x)) ;

% create vector of image points using fi function
y=[];
for j = 1:1:n
    y = [y ; fi(x(j,:),q)] ;
end
    
