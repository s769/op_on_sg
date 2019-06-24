function onbplots(T,Num)
%plots the ONB polynomials  (each OP is orthogonal to its rotations) 

%uses the functions
% SGorthoPolyspk
% SGorthoPolyssk
% gaskplot

for j=1:Num
    cla
    Q=SGorthoPolyspk(T,j-1);
    S=SGorthoPolyssk(T,j-1);
    gaskplot(sqrt(2/3)*Q + sqrt(1/3)*S,7);
    
    %formatting plot
    set(gca,'ZLim',[-10,10]);
    title(strcat('\phi_{',num2str(j-1),'}'),'fontsize',40);
    set(gca,'FontSize',20);
    
    %save each plot as an eps file
    saveas(gcf,strcat('phi',num2str(j-1)),'eps');
    %save each plot as an png file
    saveas(gcf,strcat('phi',num2str(j-1)),'png');
end
