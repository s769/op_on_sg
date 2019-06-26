function y = SG(m)
% from http://www.math.cornell.edu/~mhall/SGprograms/SG.m
% helper function for gaskplot

q0 = [cos(5*pi/12) sin(5*pi/12)]; q1=[0 0]; q2=[cos(pi/12) sin(pi/12)];
%q0 = [1/2 sqrt(3)/2]; q1 = [0, 0]; q2 = [1, 0];
y = [q0; q1; q2];

for i = 1:m
    y = [qcontract(y,q0); qcontract(y,q1); qcontract(y,q2)];
end

%plot(y(:,1),y(:,2),'.')
