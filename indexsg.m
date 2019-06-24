function n = indexsg(v,m)
% This function finds the index in TLR ordering corresponding to a
% given address

if(length(v) == 1)
% case m = 0
  n = v(1) ;
else
% general recursive formula
  n = 3^m*(v(1)-1) + indexsg(v(2:length(v)),m-1) ;
end
