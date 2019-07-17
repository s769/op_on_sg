function monomialplots(T,Num)

for j=1:Num
    clf('reset') 
    gaskplot(T(:,j,1),7);
    title(strcat('P_{',num2str(j-1),'1}'), 'FontSize', 40);
    set(gca,'FontSize',20);
    %save each plot as an jpg file
    saveas(gca,strcat('monomial',num2str(j-1),'1'),'png');
    
    clf('reset') 
    gaskplot(T(:,j,2),7);
    title(strcat('P_{',num2str(j-1),'2}'), 'FontSize', 40);
    set(gca,'FontSize',20);
    %save each plot as an jpg file
    saveas(gca,strcat('monomial',num2str(j-1),'2'),'png');
    
    clf('reset') 
    gaskplot(T(:,j,3),7);
    title(strcat('P_{',num2str(j-1),'3}'), 'FontSize', 40);
    set(gca,'FontSize',20);
    %save each plot as an jpg file
    saveas(gca,strcat('monomial',num2str(j-1),'3'),'png');
end

