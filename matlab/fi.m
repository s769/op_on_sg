function xnew = fi(x,qi)
% This function is a contractive similarity of the plane centered
% at the point qi of dilation factor 1/2.
%
% support functions: none
%
% INPUT:
% x - point in the plane
% qi - point toward which to contract distances by 1/2
%
% OUTPUT:
% evaluates the similarity
%
% from http://www.math.cornell.edu/~mhall/SGprograms/fi.m

xnew = qi + 0.5*(x-qi) ;