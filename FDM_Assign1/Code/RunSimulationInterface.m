function [u, v] = RunSimulationInterface(C1, C2, A, gridDim, deltaT, endT, x_l, x_r, MakeSBPOperators, MakeBoundariesDBC)
    rr = 0.1;
    x = linspace(x_l, x_r, gridDim)';
    [u10, u20] = InitialData(x+0.5, rr);
    
    gridDimL = ceil(gridDim/2);
    gridDimR = gridDim - gridDimL;

    u = [u10(1:gridDimL); u20(1:gridDimL)];
    v = [u10((gridDimL+1):gridDim); u20((gridDimL+1):gridDim)];

    MLeft = MakeAdvanceMatrixLeft(C1, A, gridDimL, MakeSBPOperators, MakeBoundariesDBC);
    MRight = MakeAdvanceMatrixRight(C2, A, gridDimR, MakeSBPOperators, MakeBoundariesDBC);
    [ItLeft, ItRight] = MakeInterfaceBoundaries(A, gridDimL, gridDimR, MakeSBPOperators);

    % Do iterations
    for i=1:ceil(endT/deltaT)
        ivLeft = [u(gridDimL); u(2*gridDimL)] - [v(1); v(gridDimR+1)];
        ivRight = -ivLeft;
        AdvanceLeft = @(t, u) MLeft*u + ItLeft*ivLeft;
        AdvanceRight = @(t, v) MRight*v + ItRight*ivRight;
        u = RungeKutta4Step(AdvanceLeft, i*deltaT, u, deltaT);
        v = RungeKutta4Step(AdvanceRight, i*deltaT, v, deltaT);

        uv = [u(1:gridDimL); v(1:gridDimL); u(gridDimL+1:end); v(gridDimL+1:end)];
        figure(9);
        plot(x, uv(1:gridDim), 'b', x, uv(gridDim+1:end), '--r');
        ylim([-1, 1]);
        title(['T=', num2str(i*deltaT)]);
        drawnow;
    end
end