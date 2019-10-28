close all;
clear all;


gridDim = 201;
numberOfIterations = 1000;
deltaT = 0.4/numberOfIterations;
a=1;
c=2;
%eps=.00001;
x = linspace(-1,1,gridDim)';
u = AnalyticSolution(x,0, c, a, eps);
plot(u,'linewidth',2)
title('Analytic solution: t=0')
print('ANM_Ass_3_FDM_analyticSolutionBegin','-djpeg')
for i=1:numberOfIterations
    u = AnalyticSolution(x,deltaT*i, c, a, eps);
    
    plot(u)
    ylim([0,5])
    drawnow


end %for, i
title('Analytic solution: t=T')
print('ANM_Ass_3_FDM_analyticSolutionEnd','-djpeg')


%v = RunSimulation3(gridDim, deltaT, numberOfIterations, @MakeSBP6Operators_Variable);

%plot(v);
%title('SBP 6 solution: t=T')


%% Task 4
vAnalytic = AnalyticSolution(x,0.4, c, a, eps);
v2 = RunSimulation3(gridDim, deltaT, numberOfIterations, @MakeSBP2Operators_Variable, eps);

norm(vAnalytic-v2,2)

v4 = RunSimulation3(gridDim, deltaT, numberOfIterations, @MakeSBP4Operators_Variable, eps);
norm(vAnalytic-v4,2)


v6 = RunSimulation3(gridDim, deltaT, numberOfIterations, @MakeSBP6Operators_Variable, eps);
norm(vAnalytic-v6,2)

subplot(3,1,1)
plot(v2)
title('SBP 2 solution: t=T','linewidth',1.3);
subplot(3,1,2)
plot(v4)
title('SBP 4 solution: t=T','linewidth',1.3);
subplot(3,1,3)
plot(v6)
title('SBP 6 solution: t=T','linewidth',1.3);


print('ANM_Ass_3_FDM_SBPFinalState','-djpeg');

%% Task 5
m = 101;
eps = 0.1;
v_eps1 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP4Operators_Variable, eps);

m = 101;
eps = 0.01;
v_eps2 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP4Operators_Variable, eps);

m = 101;
eps = 0.001;
v_eps3 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP4Operators_Variable, eps);

m = 101;
eps = 10^(-6);
v_eps4 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP4Operators_Variable, eps);

subplot(4,1,1)
plot(v_eps1)
title('SBP 4 solution: eps = 0.1','linewidth',1.3);
subplot(4,1,2)
plot(v_eps2)
title('SBP 4 solution: eps = 0.01','linewidth',1.3);
subplot(4,1,3)
plot(v_eps3)
title('SBP 4 solution: eps = 0.001','linewidth',1.3);
subplot(4,1,4)
plot(v_eps4)
title('SBP 4 solution: eps = 10^{-6}','linewidth',1.3);

print('ANM_Ass_3_FDM_dependenceEps','-djpeg');


%% Task 6
m = 101;
eps = 10^(-6);
v11 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP2Operators_Variable, eps);
v21 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP4Operators_Variable, eps);
v31 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP6Operators_Variable, eps);


v12 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP2Operators_Variable, eps);
v22 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP4Operators_Variable, eps);
v32 = RunSimulation3(m, deltaT, numberOfIterations, @MakeSBP6Operators_Variable, eps);

subplot(3,2,1)
plot(v11)
title('SBP 2 sol.','linewidth',1.3);
subplot(3,2,3)
plot(v21)
title('SBP 4 sol.','linewidth',1.3);
subplot(3,2,5)
plot(v31)
title('SBP 6 sol.','linewidth',1.3);
subplot(3,2,2)
plot(v12)
title('SBP 2 sol.: Artificial dissipation','linewidth',1.3)
subplot(3,2,4)
plot(v22)
title('SBP 4 sol.: Artificial dissipation','linewidth',1.3)
subplot(3,2,6)
plot(v32)
title('SBP 6 sol: Artificial dissipation','linewidth',1.3)
print('ANM_Ass_3_FDM_AD','-djpeg')


