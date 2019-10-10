close all;
%Problem2TestGLS();
Problem2TestRV();
return;
close all;
%Problem22();
close all;
%Problem23();
close all;

function Problem2TestGLS()
    h = 1/16;

    [U, M, p, t] = RunSimulationGLS(h, @InitialBoobyFunction);
    figure;
    pdesurf(p,t,U(:,size(U, 2)))
    print(['../Plots/problem_2_1_result_gls'],'-djpeg')

    animFileName = ['../Plots/animation_gls_smooth_h'];
    MakeAnimation(p, t, U, animFileName, 1/20);
end

function Problem2TestRV()
    h = 1/16;

    [U, M, p, t] = RunSimulationRV(h, @InitialBoobyFunction);
    figure;
    pdesurf(p,t,U(:,size(U, 2)))
    print(['../Plots/problem_2_1_result_rv'],'-djpeg')

    animFileName = ['../Plots/animation_rv_smooth_h'];
    MakeAnimation(p, t, U, animFileName, 1/20);
end


function Problem22()
    h = [1/4, 1/8, 1/16, 1/32];

    % Gather L2 errors for different mesh resolutions.
    L2Errors = zeros(1, length(h));
    for i=1:length(h)
        [U, M, p, t] = RunSimulationGLS(h(i), @InitialBoobyFunction);
        nodeError = U(:,1)-U(:,end);
        L2Errors(i) = sqrt(nodeError'*M*nodeError);
    end
    
    % Plot h vs L2Error
    figure;
    subplot(1,2,1)
    loglog(h, L2Errors);
    xlabel('h')
    ylabel('L^2-error')
    % ... plot h_max vs h_max^2 on the same plot
    hold on
    loglog(h, h.^1, 'b', 'linewidth', 2)
    legend('Error function', 'h_{max} versus h_{max}^{1}')
    xlabel('h')
    ylabel('L^2-error')
    % ... plot interpolation
    coeffs = polyfit(h, L2Errors, 1);
    x = linspace(0,1/4,26);
    y = coeffs(1)*x + coeffs(2);
    subplot(1,2,2)
    plot(h, L2Errors, 'r', x, y, 'b', 'linewidth', 2)
    legend('Error function', 'Linear interpolation.')
    xlabel('h')
    ylabel('L^2-error')
    print('../Plots/problem_1_2_convergence_fit','-djpeg')
end


function Problem23()
    h = [1/8, 1/16];

    % Make animated plots.
    for i=1:length(h)
        [U, M, p, t] = RunSimulationGLS(h(i), @InitialCylinderFunction);
        figure;
        pdesurf(p,t,U(:,size(U, 2)))
        print(['../Plots/problem_1_3_result_', num2str(i)],'-djpeg')

        animFileName = ['../Plots/animation_shock_h', num2str(i)];
        MakeAnimation(p, t, U, animFileName, 1/20);
    end

    % Gather L2 errors for different mesh resolutions.
    h = [1/4, 1/8, 1/16, 1/32];
    L2Errors = zeros(1, length(h));
    for i=1:length(h)
        [U, M, p, t] = RunSimulationGLS(h(i), @InitialCylinderFunction);
        nodeError = U(:,1)-U(:,end);
        L2Errors(i) = sqrt(nodeError'*M*nodeError);
    end

    % Plot h vs L2Error
    figure;
    subplot(1,2,1)
    loglog(h, L2Errors);
    xlabel('h')
    ylabel('L^2-error')
    % ... plot h_max vs h_max^2 on the same plot
    hold on
    loglog(h, h.^1, 'b', 'linewidth', 2)
    legend('Error function', 'h_{max} versus h_{max}^{1}')
    xlabel('h')
    ylabel('L^2-error')    
    % ... plot interpolation    
    coeffs = polyfit(h, L2Errors, 1);
    x = linspace(0,1/4,26);
    y = coeffs(1)*x + coeffs(2);
    subplot(1,2,2)
    plot(h, L2Errors, 'r', x, y, 'b', 'linewidth', 2)
    legend('Error function', 'Linear interpolation.')
    xlabel('h')
    ylabel('L^2-error')    
    print('../Plots/problem_1_3_errors','-djpeg')
end


function [U, M, p, t] = RunSimulationGLS(meshSize, InitialData)
    endTime = 1;
    CFL = 0.5;
    fPrimeMax = 2*pi*0.5; % ||2*pi*[-y, x]|| over a circle w/ r=0.5
    timeStep = CFL*meshSize / fPrimeMax;
    numIters = ceil(endTime/timeStep);
    timeStep = endTime / numIters;

    geometry = @circleg;
    [p,e,t] = initmesh(geometry, 'hmax', meshSize);

    xi = CreateInitialData(p, InitialData);

    delta = 0.03;
    M = MassMatrixGFEM(p,t) + delta*MassMatrixLS(p,t);
    C = ConvectionMatrixGFEM(p,t) + delta*ConvectionMatrixLS(p,t);

    U = SolverCN(M, C, xi, timeStep, numIters);
end


function [U, M, p, t] = RunSimulationRV(meshSize, InitialData)
    endTime = 1;
    CFL = 0.5;
    fPrimeMax = 2*pi*0.5; % ||2*pi*[-y, x]|| over a circle w/ r=0.5
    timeStep = CFL*meshSize / fPrimeMax;
    numIters = ceil(endTime/timeStep);
    timeStep = endTime / numIters;

    geometry = @circleg;
    [p,e,t] = initmesh(geometry, 'hmax', meshSize);

    xi = CreateInitialData(p, InitialData);

    [U, M] = SolverRV(p, t, xi, timeStep, numIters); % TODO: add proper numIters
end
