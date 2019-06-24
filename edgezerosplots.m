function edgezerosplots(T,Num,flag)
%Plots the zeros along the edges of the OP_j for j=0,...,Num-1
% if flag = 1, then using antisymm OP
% if flag = 2, then using symm OP
%
% Calls on the functions:
% edgezeros4
clf('reset')
for j=1:2:Num
    %figure
    clf('reset')
    edgezeros4(T,j,flag);
    if flag == 1
        filestring='images/p';
    else
        filestring='images/s';
    end
    %save plot as eps file
    saveas(gcf,strcat(filestring,num2str(j-1),'edgezeros'),'png');
end