function interlacing(deg,level,T)

bottomedgeindices = SGedge23(level);
load sob20coefs ops
W1 = ops;

load leg20coefs ops
W2 = ops;

% load sobolevsymmetric20 sobsymmtricops20
% W2 = sobsymmtricops20;

m=level;
%Generating legendre values
%coefflegendre=W2(deg,:);
coefflegendre=W2(deg+1,:);
legendrevalues=zeros(3^(m+1),1);
%legendrevalues = SGorthoPolyspk(T,deg);
for k=0:deg
    
        legendrevalues = legendrevalues + coefflegendre(k+1)*T(:,k+1,3);
    
end

%Generating sobolev values 
coeffsobolev=W1(deg,:);

sobolevvalues=zeros(3^(m+1),1);

for k=0:deg
    
        sobolevvalues = sobolevvalues + coeffsobolev(k+1)*T(:,k+1,3);
    
end

%evaluating legendre on edge

edgelegendre = zeros(length(bottomedgeindices),1);
in=1;
for j=1:3^(m+1)
    if j == bottomedgeindices(in)
        edgelegendre(in) = legendrevalues(j);
        in = in+1;
    end
end

%evaluating sobolev on edge

edgesobolev = zeros(length(bottomedgeindices),1);
in=1;
for j=1:3^(m+1)
    if j == bottomedgeindices(in)
        edgesobolev(in) = sobolevvalues(j);
        in = in+1;
    end
end
%clf('reset')

%Now for our next trick, we have...zeros! All the code has been copied.
%Something original. 

%z = zeros(level+1);
N = 2^(level);
w = -ones(2,2*(2^(level) + 1)-1);

for k=1:2
    

   if k == 1
        edge=edgelegendre;
        legstring = 'P_{';
       
    else
        edge=edgesobolev;
        sobstring = 'S_{';
        
    end
    %find zeros by evaluating where the OP changes sign along the boundary
    i=1;
    for j=1:N
        if edge(j)*edge(j+1) < 0 || edge(j) == 0
            %z(k,i) = (j+.5)/N;
            w(k,2*j)=2*j;
            
        end
    end
end
w(w==-1) = nan;
x=zeros(2*(2^(level) + 1)-1,1);



%formatting the plot

% set(gca,'Units','pixels');
% set(gcf,'Units','pixels');
% set(gcf, 'PaperPosition', [0 0 6 1.25]);
% set(gca,'Position',[50,25,400,75]);

figure; 
z=zeros(size(bottomedgeindices));
% plot(real(10^(3)*edgesobolev), 'r-','LineWidth',1,'DisplayName','Sobolev');
% hold on
% plot(real(edgelegendre), 'b-','LineWidth',1,'DisplayName','Legendre');

plot(1000*real(edgesobolev), 'r-','LineWidth',1,'DisplayName','Sobolev');
hold on
plot(real(edgelegendre), 'b-','LineWidth',1,'DisplayName','Legendre');


plot(z,'LineWidth',1,'Color','Black','Handlevisibility','off');
set(gca, 'XLim',[0,max(size(bottomedgeindices))-1]);

ax1 = gca; 
set(ax1,'FontSize',15);
set(ax1,'TickLabelInterpreter','latex');
legend('Location','SouthWest','Interpreter','latex');
title(strcat("Degree = ",num2str(deg-1)),'Interpreter','latex','FontSize',40); 
xlabel('Bottom Edge', 'Interpreter','latex');
ax2 = axes('Position',[0.25 0.5 0.5 0.2]);


box on;
    plot(w(1,:),x,'b*','DisplayName',strcat(legstring,num2str(deg),'}'));
    hold on
    plot(w(2,:),x,'ro','DisplayName',strcat(sobstring,num2str(deg),'}'));
    plot(1:max(size(x)),x,'k-','HandleVisibility','off');
    title('Zeroes','Interpreter','latex');
    set(gca,'XLim',[0,max(size(x))]);
    set(gca,'YLim',[-0.01,0.01]);
    set(gca,'FontSize',15);
    set(ax2,'TickLabelInterpreter','latex');
    %legend({strcat(legstring,num2str(deg),'}'),strcat(sobstring,num2str(deg),'}'),'zero'},'Location','EastOutside');
    %legend();

    set(ax2,'XTick',[]);

%saveas(gca,strcat('sd',num2str(deg-1)),'png');
    %set(gca,'YTick',-1:4:3);

