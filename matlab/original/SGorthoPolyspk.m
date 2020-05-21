function q = SGorthoPolyspk(T,deg)
%This function calculates the values of the antisymmetric OP (normalized)
%calls on function Omegas3, which imports the mathematica file with the
%values of the coefficents of the polynomials Pj

%returns vector with values of OP_level up to order m=7

m=7;
W=getOmegas3(deg);
coeff=W(deg,:);

q=zeros(3^(m+1),1);

for k=0:deg
    
        q = q + coeff(k+1)*T(:,k+1,2);
    
end

