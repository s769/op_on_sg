function w = alternateaddress(v,m)
% This function finds the alternative address for a given
% point in an SG ordering. If the point is on the boundary of SG,
% (and hence has only one such address) the same address vector is
% returned.
%
% support functions: none
%
% INPUT: 
% v - address of the point (length = m+1)
% m - the number of levels (meaning we have 3^(m+1) points) 
%
% OUTPUT: 
% w - other address vector for the same point
%
% e.g. 
% v = [1 2 3] m = 2 returns [1 3 2]
% v = [2 1 1] m = 2 returns [1 2 2]

% initialize
w = v ;

if(m == 0) % if at level zero, nothing to be done.
  return ;
else  
  % find the last entry of w which is not equal to w's final entry
  i = w(m+1) ;
  j = w(m) ;
  k = m ;
  while(j == i && k > 1)
    k = k-1 ;
    j = w(k) ;
  end
  
  % if there is such an entry, interchange its value with that of
  % all subsequent entries
  if(j ~= 0)
    w(k) = i ;
    for d = k+1:m+1
      w(d) = j ;
    end
  end
end