function q = SGorthoPolyssk(T,level)
%This function calculates the values of the symmetric OP (normalized)
%calls on function Omegas5, which imports the mathematica file with the
%values of the coefficents of the polynomials Pj

%returns vector with values of OP_level up to order m=7

m=7;
% Omegas4 - un-normalized sk
%W=getOmegas4(level);
% Omegas5 - normalized sk
W=getOmegas5(level);
coeff=W(level+1,:);

q=zeros(3^(m+1),1);

for k=0:level
    if coeff(k+1) ~= 0
        q = q + coeff(k+1)*(T(:,k+1,1)+SGrotate1(T(:,k+1,1))++SGrotate2(T(:,k+1,1)));
    end
end

