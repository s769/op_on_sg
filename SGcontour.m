function y = SGcontour(u,m)
% from http://www.math.cornell.edu/~mhall/SGprograms/SG.m
%
% plots a countour plot of the function u
%

q0 = [cos(5*pi/12) sin(5*pi/12)]; q1=[0 0]; q2=[cos(pi/12) sin(pi/12)];
y = [q0; q1; q2];

for i = 1:m
    y = [qcontract(y,q0); qcontract(y,q1); qcontract(y,q2)];
end

hold on
for j=1:3^(m+1)
    if u(j) == 0
        plot(y(j,1),y(j,2),'.','Color','black');
    elseif u(j)> 0
        plot(y(j,1),y(j,2),'.','Color','blue');%[.6 .6 .6]);
    elseif u(j) < 0
        plot(y(j,1),y(j,2),'.','Color','cyan');
    end
    if zeroflag(u,j,m) == 1
        plot(y(j,1),y(j,2),'.','Color','black');
    end
end

