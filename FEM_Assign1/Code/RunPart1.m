Problem11();
Problem12();
Problem13();
%%

function Problem11()
    h = [1/8, 1/16];

    for i=1:length(h)
        [U, M, p, t] = RunSimulationGFEM(h(i), @InitialBoobyFunction);
        figure;
        pdesurf(p,t,U(:,size(U, 2)))
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
    loglog(h, L2Errors);
    xlabel('h')
    ylabel('L^2-error')
    print('../Plots/problem_1_2_errors','-djpeg')

    % Plot parabola    
    coeffs = polyfit(h, L2Errors, 2);
    x = linspace(0,1/4,26);
    y = coeffs(1)*x.^2 + coeffs(2)*x + coeffs(3);

    figure;
    plot(h, L2Errors, 'r', x, y, 'b', 'linewidth', 2)
    legend('Error function', 'Quadratic interpolation.')
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
    loglog(h, L2Errors);
    xlabel('h')
    ylabel('L^2-error')
    print('../Plots/problem_1_3_errors','-djpeg')
end


function [U, M, p, t] = RunSimulationGFEM(meshSize, InitialData)
    endTime = 1;
    CFL = 0.5;
    fPrimeMax = 2*pi*0.5; % 2*pi*[-y, x] over a circle w/ r=0.5
    timeStep = CFL*meshSize / fPrimeMax;
    numIters = ceil(endTime/timeStep);

    geometry = @circleg;
    [p,e,t] = initmesh(geometry, 'hmax', meshSize);

    xi = CreateInitialData(p, InitialData);

    M = MassMatrixGFEM(p,t);
    C = ConvectionMatrixGFEM(p,t);

    U = SolverCN(M, C, xi, timeStep, numIters);
end
