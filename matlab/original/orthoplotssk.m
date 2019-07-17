function orthoplotssk(T,Num)
%plots the antisymmetric orthogonal polynomials 

%uses the functions
% SGorthoPolyssk
% gaskplot

for j=1:Num
%     
%     
    figure
    clf('reset')
    p = SGorthoPolyssk(T,j);
    gaskplot(p,7);
    
    %formating plot
    %set(gca,'ZLim',[-10,10]);
    title(strcat('D_{',num2str(j-1),'}'), 'FontSize', 40);
    set(gca,'FontSize',20);
    
%     % save each plot as an eps file 
%     saveas(gca,strcat('s',num2str(j-1)),'eps');
%     % save each plot as a jpg file 
     saveas(gca,strcat('Ss',num2str(j-1)),'jpg');
end

