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

[H, Q, D1, D2, M] = MakeSBP6Operators(gridDim);
[Aplus, Aminus] = SplitMatrix(A);
 
e_1=zeros(gridDim,1);e_1(1)=1;
e_m=zeros(gridDim,1);e_m(gridDim)=1;


% Prepare simulation
[u10, u20] = InitialData(x, rr);
v = [u10; u20];

Cmod = kron(C, eye(gridDim));
M = kron(A, D1);
SAT_R_CBC = -kron(Aplus,inv(H)*e_m)*kron(eye(2),e_m');
SAT_L_CBC = kron(Aminus,inv(H)*e_1)*kron(eye(2),e_1');

tau_l = [0; 1];
tau_r = [0; -1];
SAT_L_DBC = kron(tau_l, inv(H)*e_1)*kron([1 0], e_1');
SAT_R_DBC = kron(tau_r, inv(H)*e_m)*kron([1 0], e_m');

AdvanceMatrix = inv(Cmod)*(M + SAT_L_DBC + SAT_R_DBC);
eigenValues = eig(AdvanceMatrix);
fprintf("Largest real part of eig(A): %g\n", max(real(eigenValues)));
%return;

Advance = @(t, v) AdvanceMatrix*v;


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



function Am = MakeAdvanceMatrix(C, A, gridDim, MakeSBPOperators, MakeBoundaries)
    [H, Q, D1, D2, M] = MakeSBP6Operators(gridDim);

    e_1=zeros(gridDim,1);e_1(1)=1;
    e_m=zeros(gridDim,1);e_m(gridDim)=1;
end