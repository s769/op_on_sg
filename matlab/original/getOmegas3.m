function W = getOmegas3(level)
%This function  reads in the values of the coefficents of the Pj for the
%normalized antisymmetric OP (level indicates the number of OP)
% flag indicates whether to take sobolev or legendre
load sob64coefsfromNautilus ops
W = ops;

% load leg20coefs ops
% W2 = ops;
% W = [];
% for i = 2:2:size(W1,1)*2
%     W = [W;W1(i/2,:)];
%     W = [W;W2(i/2,:)];
% end
% 
% W = [W, zeros(size(W))];
    
W=W(1:level+1,1:level+1);