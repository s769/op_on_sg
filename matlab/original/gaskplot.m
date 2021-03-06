function gaskplot(f,m)
% This function plots a function defined on the level m vertices of
% SG.
%
% support functions: SG.m
%
% INPUT:
% f - vector of size 3^(m+1) holding function values on vertices
% m - level of the function
%
% OUTPUT:
% plots figure of the graph
%
% from http://www.math.cornell.edu/~mhall/SGprograms/gaskplot.m

% get coordinates of midpoints of edges from SGedges function
y = SG(m) ;

%hold off

% the variable s allows for the flexibility of handling either a
% row or column vector. we check which case it is using size(f) and
% do the appropriate plot
s = size(f) ;


    
    if s(2)==1
      for k=1:1:3^m

        plot3([y(3*k-2:3*k,1); y(3*k-2,1)],...
          [y(3*k-2:3*k,2) ;y(3*k-2,2)],...
          [f(3*k-2:3*k);  f(3*k-2)])
        hold on
      end
    else
      for k=1:1:3^m

        plot3([y(3*k-2:3*k,1); y(3*k-2,1)],...
          [y(3*k-2:3*k,2); y(3*k-2,2)],...
          [f(3*k-2:3*k)  f(3*k-2)])
        hold on
      end
    end
    

    