% Model definition: Cu_t  = Au_x
permittivity = 1;
permeability = 1;

A = [0 1; 1 0];
C = [permittivity 0; 0 permeability];
     
     
% Simulation definition     
gridDim = 129;
x_l = -1;
x_r = 1;
x = linspace(x_l, x_r, gridDim)';
rr = 0.1;

[H, Q, D1, D2, M] = MakeSBPOperators(gridDim);
[Aplus, Aminus] = SplitMatrix(A);
 
e_1=zeros(gridDim,1);e_1(1)=1;
e_m=zeros(gridDim,1);e_m(gridDim)=1;


% Prepare simulation
[u10, u20] = InitialData(x, rr);
v = [u10; u20];

Cmod = kron(C, eye(gridDim));
P = kron(A, D1);
SATs = 0;
Advance = @(t, v) Cmod\(P*v + SATs);


% Do iterations
h = 1/60;
for i=1:60
    v = RungeKutta4Step(Advance, i*h, v, h);
    u1 = v(1:gridDim);
    u2 = v(gridDim+1:end);
    plot(x, u1, 'r', x, u2, '--b', 'LineWidth', 2);
    ylim([-1, 1]);
    drawnow;
    pause(0.1);
end
