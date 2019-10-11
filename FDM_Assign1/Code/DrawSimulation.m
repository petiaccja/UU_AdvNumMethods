function DrawSimulation(C, A, gridDim, deltaT, endT, x_l, x_r)
    rr = 0.1;
    x = linspace(x_l, x_r, gridDim)';
    domainLength = x_r - x_l;

    % Prepare simulation
    [u10, u20] = InitialData(x, rr);
    v = [u10; u20];

    M = MakeAdvanceMatrix(C, A, gridDim, @MakeSBP6Operators, @MakeBoundariesCBC);
    eigenValues = eig(M);
    fprintf('Largest real part of eig(A): %g\n', max(real(eigenValues)));

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