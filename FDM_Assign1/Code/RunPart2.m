close all;
clear all;

% Model definition: Cu_t  = Au_x
permittivity1 = 0.3;
permittivity2 = 1;
permeability = 1;

A = [0 1; 1 0];
C1 = [permittivity1 0; 0 permeability];
C2 = [permittivity2 0; 0 permeability];


gridDim = 202;

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

refractiveIndex1 = sqrt(permittivity1);
refractiveIndex2 = sqrt(permittivity2);

Tanalytic = abs(2*refractiveIndex1/(refractiveIndex1+refractiveIndex2))
Ranalytic = abs((refractiveIndex1-refractiveIndex2)/(refractiveIndex1+refractiveIndex2))

c1 = 1/refractiveIndex1;
c2 = 1/refractiveIndex2;
xWaveOriginal = -1 + (c1*endT - 0.5);
xWaveReflected = 0 - (c1*endT - 0.5);
xWaveTransmitted = 0 + (endT - 0.5/c1)*c2;

% analytic solution for t=0.42
rr = 0.1;
x = linspace(x_l, x_r, gridDim)';
[u1o, u2o] = InitialData(x-xWaveOriginal, rr);
[u1t, u2t] = InitialData(x-xWaveTransmitted, rr*c2/c1);
[u1r, u2r] = InitialData(x-xWaveReflected, rr);
u1end = 0.5*u2o + 0.5*Ranalytic*u2r + 0.5*Tanalytic*u2t/refractiveIndex1;
u2end = -0.5*u2o/refractiveIndex1 + 0.5*Ranalytic*u2r/refractiveIndex1 + -0.5*Tanalytic*u2t/refractiveIndex1;

figure;
plot(x, u1end, x, u2end);


[vl, vr] = RunSimulationInterface(C1, C2, A, gridDim, deltaT, endT, x_l, x_r, @MakeSBP6Operators, @MakeBoundariesDBC);
gridDimL = length(vl)/2;
gridDimR = length(vr)/2;
v = [vl(1:gridDimL); vr(1:gridDimL); vl(gridDimL+1:end); vr(gridDimL+1:end)];
figure;
plot(x, v(1:gridDim), 'b', x, v(gridDim+1:end), '--r');
%ylim([-1, 1]);
xlabel('x');
ylabel('u(x, t)');
legend('E=u_1', 'H=u_2');
title('Solution: Dirichlet BCs, 6th order, m=201');
print('FDM_Ass1_DBC_6thOrder','-djpeg')



amplitudeOriginal = -min(vl(1:gridDimL));
amplitudeReflected = max(vl(1:gridDimL));
amplitudeTransmitted = -min(vr(1:gridDimR));

Tnumeric = amplitudeTransmitted/amplitudeOriginal
Rnumeric = amplitudeReflected/amplitudeOriginal









