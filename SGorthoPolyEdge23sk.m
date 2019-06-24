function edge = SGorthoPolyEdge23sk(T,level)
%
% evaluates and plots the orthonormal polynomial P_level_base along the
% edge between v1 and v3
%
% calls on the functions:
% SGedge23
% SGorthoPolyssk

m=7;
indices = SGedge23(m);
Q = SGorthoPolyssk(T,level);

edge = zeros(2^(m+1),1);
in=1;
for j=1:3^(m+1)
    if j == indices(in)
        edge(in) = Q(j);
        in = in+1;
    end
end

