function l = alternateindex(k,m)
% This function finds the alternative index for a given
% point in TLR ordering on SG. If the point is on the boundary of SG,
% (and hence has only one such index) the same address vector is
% returned.
%
% support functions: alternateaddress.m
%
% INPUT: 
% k - index of the point
% m - the number of levels (meaning we have 3^(m+1) possible indices) 
%
% OUTPUT: 
% l - other index for the same point
%
% e.g. 
% alternateindex(2,1) returns 5
% alternateindex(3,1) returns 7


% We convert to an address, go through the alternate address
% function, and convert back to an index
v = address(k,m) ;
w = alternateaddress(v,m) ;
l = indexsg(w,m) ;