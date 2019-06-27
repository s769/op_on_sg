function trialplots(n)

x = linspace(0,2*pi,100)

for i = 1:n
  
subplot(n,1,i)
plot(x, sin(i*x))
hold on 
plot(x,zeros(size(x)))

%formatting plot
    %set(gca,'XLim',[0,length(edge)]);
    %set(gca,'XTick', []);
    xlabel('xlabel','FontSize',30);
    %set(gca,'YLim',[-4,4]);
    %title(strcat(legstring,num2str(j-1),'}'), 'FontSize',40);
    %set(gca,'FontSize',30);
end

