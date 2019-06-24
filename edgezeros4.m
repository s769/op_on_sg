function z = edgezeros4(T,level1,flag)
%Plots the zeros on the bottom edge of the OP_level1
% if flag = 1, then looking at antisymm OP
% if flag = 2, then looking at symm OP
%
% Calls on the functions:
% SGorthoPolyEdge23pk
% SGorthoPolyEdge23sk

z = zeros(level1+1);
N=256;
w = -ones(2,N);

for k=1:2
    if flag ==1
        edge=SGorthoPolyEdge23pk(T,level1+k-1);
        legstring = 'p_{';
    else
        edge=SGorthoPolyEdge23sk(T,level1+k-1);
        legstring = 's_{';
    end
    
    %find zeros by evaluating where the OP changes sign along the boundary
    i=1;
    for j=1:N-1
        if edge(j)*edge(j+1) < 0 || edge(j) == 0
            z(k,i) = (j+.5)/N;
            w(k,j)=1;
            i=i+1;
        end
    end
end


%formatting the plot
x=1:N;
set(gca,'Units','pixels');
set(gcf,'Units','pixels');
set(gcf, 'PaperPosition', [0 0 6 1.25]);
set(gca,'Position',[50,25,400,75]);
plot(x,w(1,:),'o',x,w(2,:),'*');
set(gca,'XLim',[0,N]);
%set(gca,'YLim',[0,2]);
legend(strcat(legstring,num2str(level1-1),'}'),strcat(legstring,num2str(level1),'}'),'Location','EastOutside');
legend('boxoff');
%title('Location of zeros for each Qk on the bottom edge');
set(gca,'XTick',-N:3*N:2*N);
%set(gca,'YTick',-1:4:3);
