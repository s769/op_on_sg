function edgeplots(T,Num,flag)

%
%plots the polynomials along the edges, plots them all together
% flag = 1 plots the edge of the antisymm OP along the bottom edge
% flag = 2 plots the edge of the symm OP along the bottom edge
% flag = 3 plots the edge of the antisymm OP along the side edge
clf('reset')

for j=1:Num
    %cla
    if flag==1
        edge=SGorthoPolyEdge23pk(T,j);
        titlestring = 'bottomedge';
        xlabelstr = 'bottom edge';
        legstring = 'S_{';
        filestring = 'images/s';
    elseif flag==2
        edge=SGorthoPolyEdge23sk(T,j);
        titlestring = 'bottomedge';
        xlabelstr = 'bottom edge';
        legstring = 'D_{';
        filestring = 'images/s';
    elseif flag==3
        edge=SGorthoPolyEdge13pk(T,j);
        titlestring = 'sideedge';
        xlabelstr = 'side edge';
        legstring = 'S_{';
        filestring = 'images/s';
    elseif flag==4
        edge=SGorthoPolyEdge23phik(T,j);
        titlestring = 'bottomedge';
        xlabelstr = 'bottom edge';
        legstring = '\phi_{';
        filestring = 'images/phi';
    else
        edge=SGorthoPolyEdge13phik(T,j);
        titlestring = 'sideedge';
        xlabelstr = 'side edge';
        legstring = '\phi_{';
        filestring = 'images/phi';        
    end
    figure
    
    %subplot(Num,2,j)
    clf('reset')
    z=zeros(size(edge));
    plot(real((edge)), 'LineWidth',1);
    hold on
    plot(z,'LineWidth',1,'Color','Black');
    
    
    %formatting plot
    set(gca,'XLim',[0,length(edge)]);
    set(gca,'XTick', []);
    xlabel(xlabelstr,'FontSize',30);
    %set(gca,'YLim',[-4,4]);
    title(strcat(legstring,num2str(j-1),'}'), 'FontSize',40);
    set(gca,'FontSize',30);
    
    %save each plot as an eps file
    %saveas(gca,strcat(filestring,num2str(j-1),titlestring),'eps');
    %save each plot as an jpg file
    saveas(gca,strcat(filestring,num2str(j-1),titlestring),'jpg');
end

