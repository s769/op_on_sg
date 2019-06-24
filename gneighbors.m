function nbhd = gneighbors(v,m)
% This function finds the addresses of those edges which are
% incident to a given vertex. If the given vertex is on the
% boundary of SG we return two copies of each edge. In fact the
% addresses are the same as those of the adjacent vertices (compare
% vneighbors.m), but we correct the ordering here. As a result, the
% returned matrix has as rows those addresses corresponding to the
% edges between this vertex and the vertex whose address is found
% in the corresponding row of the matrix returned by vneighbors.m
%
% support functions: alternateaddress.m
%
% INPUT:
% v - address of a level m vertex
% m - the level of SG graph
%
% OUTPUT:
% nbhd - matrix with 4 rows giving the addresses of edges touching
% this vertex
%
% e.g.
%      /\
%     /__\
%    /\  /\
%   /__\/__\
%
% gneighbors([1 2]) = [1 1; 1 3; 2 3; 2 2]
%


% find alternate address of our vertex
w = alternateaddress(v,m) ;

one = ones(size(v)) ;
di = zeros(size(v)) ;
di(length(di)) = 1 ;

% you can check the following formulas are correct: given the
% address of a point corresponding to that point's position within
% a given cell, the incident edges in that cell have the two
% other vertex addresses corresponding to that cell
nbhd(1,:) = mod(v-one+2*di,3)+one ;
nbhd(2,:) = mod(v-one+1*di,3)+one ;
nbhd(3,:) = mod(w-one+2*di,3)+one ;
nbhd(4,:) = mod(w-one+1*di,3)+one ;