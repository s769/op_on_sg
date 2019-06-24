function flag=zeroflag(u,j,m)

nbhd = gneighbors(address(j,m),m);
pflags = zeros(4,1);
nflags = zeros(4,1);
nbr = zeros(4,1);
for j=1:4
    nbr(j) = indexsg(nbhd(j,:),m);
    if u(nbr(j))>0
        pflags(j) =1;
    elseif u(nbr(j))<0
        nflags(j) = 1;
    end        
end

if pflags == not(nflags)
    if pflags == ones(4,1)
        flag=0;
    elseif pflags == zeros(4,1)
        flag=0;
    else
        flag=1;
    end
else
    flag=0;
end
