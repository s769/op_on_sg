function SGplotsOverview(T)

Num=20;

%This reads in the values of the polynomials Pj for j=0,...,19
%T is a large matrix with the values of the polynomials up to level m=7
% (This takes a while to read, if doing multiple runs of this program, might
% want to read in T and pass it to this program)
%T=readPolys();

%To plot the antisymmetric OP, do the following:
%Num is the number of OP plotted
%orthoplotspk(T,Num);
%figure

% %To plot the symmetric OP, do the following:
% %Num is the number of OP plotted
orthoplotssk(T,Num);
figure
% 
% %To plot the ONP, do the following:
% %Num is the number of OP plotted
% onbplots(T,Num);
% figure
% 
% 
% %To plot the edges of the OP, do the following:
% %Num is the number of OP plotted
% %for the bottom edge of the antisymm OP
%edgeplots(T,Num,1);
%figure
% %for the bottom edge of the symm OP
%edgeplots(T,Num,2);
% figure
% %for the side edge of the antisymm OP
edgeplots(T,Num,3)
figure
% 
% %To plot the zeros along the edges of the OP, do the following:
% %Num is the number of OP plotted
% %for the bottom edge of the antisymm OP
 %edgezerosplots(T,Num,1);
% figure
% %for the bottom edge of the symm OP
% edgezerosplots(T,Num,2);