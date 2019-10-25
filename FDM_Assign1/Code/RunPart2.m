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

% Calculate convergence rate
q2 = CalculateConvergenceInterface(C1, C2, A, @MakeSBP2Operators, @MakeBoundariesDBC)
q4 = CalculateConvergenceInterface(C1, C2, A, @MakeSBP4Operators, @MakeBoundariesDBC)
q6 = CalculateConvergenceInterface(C1, C2, A, @MakeSBP6Operators, @MakeBoundariesDBC)


% Plot results for SBP6
deltaT = 0.1/gridDim;
numIters = ceil(0.42/deltaT);
endT = deltaT*numIters; %ceil(1.75/deltaT)*deltaT;
x_l = -1;
x_r = 1;
x = linspace(x_l, x_r, gridDim)';

refractiveIndex1 = sqrt(permittivity1);
refractiveIndex2 = sqrt(permittivity2);

Tanalytic = abs(2*refractiveIndex1/(refractiveIndex1+refractiveIndex2))
Ranalytic = abs((refractiveIndex1-refractiveIndex2)/(refractiveIndex1+refractiveIndex2))


[vl, vr] = RunSimulationInterface(C1, C2, A, gridDim, deltaT, numIters, x_l, x_r, @MakeSBP6Operators, @MakeBoundariesDBC);
gridDimL = length(vl)/2;
gridDimR = length(vr)/2;
v = [vl(1:gridDimL); vr(2:gridDimR); vl(gridDimL+1:end); vr(gridDimR+2:end)];
figure;
plot(x, v(1:gridDim), 'b', x, v(gridDim+1:end), '--r');
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









