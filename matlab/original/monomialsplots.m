function monomialsplots(Num,k,level,T)
% Num = the degree up to which you want the monomials plotted. 
%Default is Num = 20. 
% k = family of polynomial (constant symmetric, symmetric, and
% antisymmetric)
% level = level of gasket. Default value should be 7 

%Functions required:
% gaskplot.m


for i = 1:Num 
    figure
    gaskplot(5.7442*T(:,i,k),level)
    %saveas(gca,strcat('monomial',num2str(i-1),'_',num2str(k)),'png');
end

