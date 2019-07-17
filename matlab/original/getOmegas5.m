function W = getOmegas5(level)
%This function  reads in the values of the coefficents of the Pj for the
%normalized symmetric OP (level indicates the number of OP)

%W=dlmread('sk-recursionN');

load sobolevsymmetric20 sobsymmtricops20
W = sobsymmtricops20;
W=W(1:level+1,1:level+1);