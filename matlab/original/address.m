function v = address(k,m)
% This function returns a vector listing the address of a 
% given index k in a level m SG-type ordering (TLR ordering).
%
% support functions: none
%
% INPUT: 
% k - the index whose address to find; 
% m - the number of levels (meaning we have 3^(m+1) points) 
%
% OUTPUT: vector of length m+1, listing off a sequence of 0's, 1's,
% and 2's, which locate the given index.
%
% e.g.
% k = 4, m = 1 returns [2 1] (second cell of three, first position)
% k = 21 m = 2 returns [3 1 3]

% initialize
v = [] ;

% compute cell positions at all levels > 0; starting with
% the highest level
while(m > 0)
  r = floor((k-1) / 3^m) ;  
  v = [v r+1] ;
  k = k - r*3^m ;
  m = m-1 ;
end

% add in k mod 3 at the end
v = [v k] ; 