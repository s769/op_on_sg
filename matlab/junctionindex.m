function v=junctionindex(add0, m)
%
% given address add0 in V_n
% this outputs the index of add0 in V_m for m>=n
%
% INPUT:
% add0 and address of the form [1 2]
% m a level greater than the length of add0
%
% OUTPUT:
% v an index of the point indicated in add0, but in Vm
%
% support functions: indexsg.m
%
% e.g.
% junctionindex([1 2], 3) = 11
%

n=length(add0)-1;
if m>=n
   k=m-n;
   add1 = add0;
   q = add1(n+1);
   for j=1:k-1
       add1 = [add1 q];
   end
else
    m=n;
    add1 = add0;
end
v=indexsg(add1,m-1);
    
