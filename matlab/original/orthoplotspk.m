function orthoplotspk(T,Num)
%plots the antisymmetric orthogonal polynomials 

%uses the functions
% SGorthoPolyspk
% gaskplot

for j=1:Num
    %if mod(j, 2) == 1
        %figure
        
    %end
    %clf('reset')
    figure('visible','off')
    p = SGorthoPolyspk(T,j);
    gaskplot(p,7);
    
    %formating plot
    %set(gca,'ZLim',auto);
%     if mod(j, 2) == 1
%         title(strcat('Q_{',num2str(floor((j-1)/2)),'}'), 'FontSize', 40);
%     end
%     if mod(j, 2) == 0
%         
%     end
    title(strcat('$T_{',num2str(j-1),'}$'), 'Interpreter','latex','FontSize', 40);
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex');
    % save each plot as an eps file 
    %saveas(gca,strcat('images/p',num2str(j-1),'.eps'),'eps');
    %save as png
    filestring = strcat('images/H1_k2_lam1_0_20/T_',num2str(j-1),'.png');
    saveas(gca,filestring,'png');
end

