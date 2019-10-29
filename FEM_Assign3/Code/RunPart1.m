Problem11();
%close all;
%Problem12();
%close all;
%Problem13();
close all;

function Problem11()
    h = [1/8];%, 1/16];

    for i=1:length(h)
        [U, M, p, t] = RunSimulationGFEM(h(i), @InitialFunction);
        figure;
        pdesurf(p,t,U(:,size(U, 2)))
        xlim([-4,2])
        ylim([-4,2])
        zlim([0,12])
        print(['../Plots/problem_1_1_result_', num2str(i)],'-djpeg')

        animFileName = ['../Plots/animation_smooth_h', num2str(i)];
        MakeAnimation(p, t, U, animFileName, 1/20);
    end
end


function Problem12()
    h = [1/4, 1/8, 1/16, 1/32];

    % Gather L2 errors for different mesh resolutions.
    L2Errors = zeros(1, length(h));
    for i=1:length(h)
        [U, M, p, t] = RunSimulationGFEM(h(i), @InitialBoobyFunction);
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


function Problem13()
    h = [1/8, 1/16];

    % Make animated plots.
    for i=1:length(h)
        [U, M, p, t] = RunSimulationGFEM(h(i), @InitialCylinderFunction);
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
        [U, M, p, t] = RunSimulationGFEM(h(i), @InitialCylinderFunction);
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

function geom = rectangleg()
    geom = [
        2 2 2   -2.5 1.5  1 0; ...
        2 2 -2  1.5 1.5   1 0; ...
        2 -2 -2 1.5 -2.5  1 0; ...
        2 -2 2  -2.5 -2.5 1 0]';
end

function [U, M, p, t] = RunSimulationGFEM(meshSize, InitialData)
    endTime = 1;
    timeStep = 0.01;
    numIters = ceil(endTime/timeStep);
    endTime = numIters*timeStep;

    geometry = @rectangleg;
    [p,e,t] = initmesh(geometry(), 'hmax', meshSize);

    xi = CreateInitialData(p, InitialData);
    figure(9);
    pdesurf(p, t, xi);
    title('Initial stuff')

    M = MassMatrixGFEM(p,t);

    U = SolverCN(M, p, t, xi, timeStep, numIters);
end
