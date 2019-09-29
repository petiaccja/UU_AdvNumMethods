close all;
clear all;

% Model definition: Cu_t  = Au_x
permittivity = 1;
permeability = 1;

A = [0 1; 1 0];
C = [permittivity 0; 0 permeability];
     
     
% Simulation definition     
gridDim = 31;
deltaT = 0.1/gridDim;
endT = ceil(0.75/deltaT)*deltaT;
x_l = -1;
x_r = 1;

q2 = CalculateConvergence(C, A, @MakeSBP2Operators, @MakeBoundariesDBC)
q4 = CalculateConvergence(C, A, @MakeSBP4Operators, @MakeBoundariesDBC)
q6 = CalculateConvergence(C, A, @MakeSBP6Operators, @MakeBoundariesDBC)

function v = RunSimulation(C, A, gridDim, deltaT, endT, x_l, x_r, MakeSBPOperators, MakeBoundariesDBC)
    rr = 0.1;
    x = linspace(x_l, x_r, gridDim)';
    [u10, u20] = InitialData(x, rr);
    v = [u10; u20];

    M = MakeAdvanceMatrix(C, A, gridDim, MakeSBPOperators, MakeBoundariesDBC);

    Advance = @(t, v) M*v;

    % Do iterations
    for i=1:ceil(endT/deltaT)
        v = RungeKutta4Step(Advance, i*deltaT, v, deltaT);
    end
end


function DrawSimulation(C, A, gridDim, deltaT, endT, x_l, x_r)
    rr = 0.1;
    x = linspace(x_l, x_r, gridDim)';
    domainLength = x_r - x_l;

    % Prepare simulation
    [u10, u20] = InitialData(x, rr);
    v = [u10; u20];

    M = MakeAdvanceMatrix(C, A, gridDim, @MakeSBP6Operators, @MakeBoundariesDBC);
    eigenValues = eig(M);
    fprintf("Largest real part of eig(A): %g\n", max(real(eigenValues)));

    Advance = @(t, v) M*v;


    % Do iterations
    figure;
    for i=1:ceil(endT/deltaT)
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
end


function q = CalculateConvergence(C, A, MakeSBPOperators, MakeBoundariesDBC)
    gridDims = [31, 61];
    x_l = -1;
    x_r = 1;
    norms = [0, 0];    
    
    for i = 1:2
        gridDim = gridDims(i);
        deltaX = 1.0/gridDim;
        deltaT = 0.1*deltaX;
        endT = ceil(0.75/deltaT)*deltaT;
        x = linspace(x_l, x_r, gridDim)';
        
        v = RunSimulation(C, A, gridDim, deltaT, endT, x_l, x_r, MakeSBPOperators, MakeBoundariesDBC);
        figure;
        plot(x, v(1:gridDim), x, v(gridDim+1:end));
        [u1, u2] = AnalyticSolution(x, endT, x_r-x_l, 0.1);
        hold on;
        plot(x, u1, x, u2);
        hold off;
        norms(i) = sqrt(deltaX)*norm(v - [u1; u2]);
    end
    
    q = log10(norms(1)/norms(2)) / log10(gridDims(1)/gridDims(2));
end
