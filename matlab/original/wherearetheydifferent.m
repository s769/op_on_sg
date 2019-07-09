function wherearetheydifferent(k)

for i = 1:20
    
    figure;
    %disp(i)
   gaskplot(S(:,i,k),7,1)
   gaskplot(Torig(:,i,k),7,0)
end

