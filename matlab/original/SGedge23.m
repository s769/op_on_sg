function indices = SGedge23(m)
%
% returns the indices of the points along the edge of the 
% Sierpenski gasket between v2 and v3 on the level m graph.
%
% calls on the functions:
% indexsg

indices = zeros(2^(m+1),1);
for j=1:2^(m+1)
    v=[];
    l=j;
    for k=1:m+1
        if mod(l,2) == 0
            v=[v 2];
        else
            v=[v 3];
        end    
        l=floor(l./2);
    end
    indices(j) = indexsg(v,m);
end
indices = sort(indices);

indices(2:2:(length(indices)-1)) = [];

