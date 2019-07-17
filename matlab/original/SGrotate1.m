function rot = SGrotate1(graph)
%
% rotates the level m graph of SG so that q1 is at the top
% used for making fully symmetric polynomial
%

N=length(graph);
m=log3(N)-1;
if m>0
    graph1 = graph(1:(N/3));
    graph2 = graph((N/3)+1:(2*N/3));
    graph3 = graph((2*N/3)+1:N);
    rot = [SGrotate1(graph2); SGrotate1(graph3); SGrotate1(graph1)];
else
    rot = [graph(2); graph(3); graph(1)];
end