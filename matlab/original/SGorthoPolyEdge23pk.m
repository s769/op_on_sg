function edge = SGorthoPolyEdge23pk(T,level)
%
% evaluates and plots the antisymmetric orthonormal polynomial
% P_level_base along the edge between v1 and v3 (bottom edge)
%
% calls on the functions:
% SGedge23
% SGorthoPolyspk

m=7;
indices = SGedge23(m);
Q = SGorthoPolyspk(T,level);

edge = zeros(length(indices),1);
in=1;
for j=1:3^(m+1)
    if j == indices(in)
        edge(in) = Q(j);
        in = in+1;
    end
end


