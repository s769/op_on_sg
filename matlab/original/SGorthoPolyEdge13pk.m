function edge = SGorthoPolyEdge13pk(T,level)
%
% evaluates and plots the orthonormal polynomial P_level_base along the
% edge between v1 and v3 (side edge)
%
% calls on the functions:
% SGedg123
% SGorthoPolyspk

m=7;
indices = SGedge13(m);
Q = SGorthoPolyspk(T,level);

edge = zeros(2^(m+1),1);
in=1;
for j=1:3^(m+1)
    if j == indices(in)
        edge(in) = Q(j);
        in = in+1;
    end
end