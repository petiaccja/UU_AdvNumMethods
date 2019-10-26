close all;
clear all;


gridDim = 101;
numberOfIterations = 1000;
deltaT = 0.4/numberOfIterations;
a=1;
c=2;
%eps=.00001;
x = linspace(-1,1,gridDim)'
% for i=1:numberOfIterations
%     u = AnalyticSolution(x,deltaT*i, c, a, eps);
%     plot(u
%     ylim([0,5])
%     drawnow
% 
% end %for, i



v = RunSimulation3(gridDim, deltaT, numberOfIterations, @MakeSBP2Operators_Variable);

v
plot(v);




