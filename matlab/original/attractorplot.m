function attractorplot(T,index,degree)

for k=1:degree-1
    kdegree = SGorthoPolyspk(T,k);
    x = kdegree(index);
    kplusonedegree = SGorthoPolyspk(T,k+1);
    y = kplusonedegree(index);
    kplustwodegree = SGorthoPolyspk(T,k+2);
    z = kplustwodegree(index);
    plot(log(y^2),log(x*z),'bo-');
    hold on
end

