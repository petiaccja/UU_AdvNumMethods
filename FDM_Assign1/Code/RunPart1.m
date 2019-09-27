close all;
clear all;

% Model definition: Cu_t  = Au_x
permittivity = 1;
permeability = 1;

A = [0 1; 1 0];
C = [permittivity 0; 0 permeability];
     
     
% Simulation definition     
gridDim = 201;
deltaX = 1/gridDim;
deltaT = 0.2*deltaX;
x_l = -1;
x_r = 1;
x = linspace(x_l, x_r, gridDim)';
domainLength = x_r - x_l;
rr = 0.1;


% Prepare simulation
[u10, u20] = InitialData(x, rr);
v = [u10; u20];

M = MakeAdvanceMatrix(C, A, gridDim, @MakeSBP6Operators, @MakeBoundariesDBC);
eigenValues = eig(M);
fprintf("Largest real part of eig(A): %g\n", max(real(eigenValues)));
%return;

Advance = @(t, v) M*v;


% Do iterations
for i=1:ceil(2/deltaT)
    v = RungeKutta4Step(Advance, i*deltaT, v, deltaT);
    u1 = v(1:gridDim);
    u2 = v(gridDim+1:end);
    [u1_an, u2_an] = AnalyticSolution(x, i*deltaT, domainLength, rr);
    plot(x, u1, 'b', x, u2, '--r', 'LineWidth', 2);
    hold on;
    plot(x, u1_an, 'g', x, u2_an, '--m', 'LineWidth', 1);
    hold off;
    ylim([-1, 1]);
    drawnow;
    pause(0.005);
end
