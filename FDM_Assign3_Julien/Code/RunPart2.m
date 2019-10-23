close all;
clear all;

% Model definition: Cu_t  = Au_x
permittivity1 = 0.3;
permittivity2 = 1;
permeability = 1;

A = [0 1; 1 0];
C1 = [permittivity1 0; 0 permeability];
C2 = [permittivity2 0; 0 permeability];


gridDim = 201;

%% Calculate convergence rate
%q2 = CalculateConvergence(C, A, @MakeSBP2Operators, @MakeBoundariesDBC)
%q4 = CalculateConvergence(C, A, @MakeSBP4Operators, @MakeBoundariesDBC)
%q6 = CalculateConvergence(C, A, @MakeSBP6Operators, @MakeBoundariesDBC)

%% Plot results for SBP6
gridDim = 201;
deltaT = 0.1/gridDim;
% gets unstable around deltaT=1.37*h for DBC w/ m=201 SBP6
% ~2.1 for SBP2
endT = .42; %ceil(1.75/deltaT)*deltaT;
x_l = -1;
x_r = 1;
x = linspace(x_l, x_r, gridDim)';

[vl, vr] = RunSimulationInterface(C1, C2, A, gridDim, deltaT, endT, x_l, x_r, @MakeSBP6Operators, @MakeBoundariesDBC);
gridDimL = length(vl)/2;
gridDimR = length(vr)/2;
v = [vl(1:gridDimL); vr(1:gridDimL); vl(gridDimL+1:end); vr(gridDimL+1:end)];
%figure;
%plot(x, v(1:gridDim), 'b', x, v(gridDim+1:end), '--r');
%ylim([-1, 1]);
%xlabel('x');
%ylabel('u(x, t)');
%legend('E=u_1', 'H=u_2');
%title('Solution: Dirichlet BCs, 6th order, m=201');
%print('FDM_Ass1_DBC_6thOrder','-djpeg')


refrectingIndex1 = sqrt(permittivity1);
refrectingIndex2 = sqrt(permittivity2);

Tanalytic = 2*refrectingIndex1/(refrectingIndex1+refrectingIndex2)
Ranalytic = (refrectingIndex1-refrectingIndex2)/(refrectingIndex1+refrectingIndex2)

amplitudeOriginal = -min(vl(1:gridDimL));
amplitudeReflected = max(vl(1:gridDimL));
amplitudeTransmitted = -min(vr(1:gridDimR));

Tnumeric = amplitudeTransmitted/amplitudeOriginal
Rnumeric = amplitudeReflected/amplitudeOriginal



