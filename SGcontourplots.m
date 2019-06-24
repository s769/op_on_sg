function SGcontourplots(T,Num,flag)

m=7;
for j=1:Num
    Q=SGorthoPolyspk(T,j-1);
    S=SGorthoPolyssk(T,j-1);
    phi=sqrt(2/3)*Q + sqrt(1/3)*S;
    hold on
        
    if flag==1
        clf('reset')
        SGcontour(Q,m);
        %formatting the graph
        axis([0 1 0 1])
        %axis square
        box
        set(gca,'XTick', []);
        set(gca,'YTick', []);
        title(strcat('contour plot of Q_{',num2str(j-1),'}'),'FontSize', 20);
        % save each plot as an eps file 
        saveas(gca,strcat('contourp',num2str(j-1),'.eps'),'epsc');
        % save each plot as an jpg file 
        saveas(gca,strcat('contourp',num2str(j-1),'.png'),'png');  
    elseif flag==2
        clf('reset')
        SGcontour(S,m);
        %formatting the graph
        axis([0 1 0 1])
        box
        set(gca,'XTick', []);
        set(gca,'YTick', []);
        title(strcat('contour plot of S_{',num2str(j-1),'}'),'FontSize', 20);
        % save each plot as an eps file 
        saveas(gca,strcat('contours',num2str(j-1),'.eps'),'epsc');
        % save each plot as an jpg file 
        saveas(gca,strcat('contours',num2str(j-1),'.png'),'png');
    else
        clf('reset')
        SGcontour(phi,m);
        %formatting the graph
        axis([0 1 0 1])
        box
        set(gca,'XTick', []);
        set(gca,'YTick', []);
        title(strcat('contour plot of \phi_{',num2str(j-1),'}'),'FontSize', 20);
        % save each plot as an eps file 
        saveas(gca,strcat('contourphi',num2str(j-1),'.eps'),'epsc');
        % save each plot as an jpg file 
        saveas(gca,strcat('contourphi',num2str(j-1),'.png'),'png');
    end
end