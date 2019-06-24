function edge = SGorthoPolyEdge23phik(T,level)
%
% evaluates and plots the orthonormal polynomial phi_level_base along the
% edge between v1 and v3 (side edge)
%
% calls on the functions:
% SGedge23
% SGorthoPolyspk
% SGorthoPolyssk

m=7;
indices = SGedge23(m);
Q = SGorthoPolyspk(T,level);
S = SGorthoPolyssk(T,level);
phi=sqrt(2/3)*Q + sqrt(1/3)*S;

edge = zeros(2^(m+1),1);
in=1;
for j=1:3^(m+1)
    if j == indices(in)
        edge(in) = phi(j);
        in = in+1;
    end
end