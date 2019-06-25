function rot = SGrotate2(graph)
%
% rotates the level m graph of SG so that q2 is at the top
% used for making fully symmetric polynomial
%

N=length(graph);
m=log3(N)-1;
if m>0
    graph1 = graph(1:(N/3));
    graph2 = graph((N/3)+1:(2*N/3));
    graph3 = graph((2*N/3)+1:N);
    rot = [SGrotate2(graph3); SGrotate2(graph1); SGrotate2(graph2)];
else
    rot = [graph(3); graph(1); graph(2)];
end