function [vl, vr] = RunSimulationInterface(C1, C2, A, gridDim, deltaT, numIters, x_l, x_r, MakeSBPOperators, MakeBoundariesDBC)
    rr = 0.1;
    x = linspace(x_l, x_r, gridDim)';
    [u10, u20] = InitialData(x+0.5, rr);
    
    gridDimL = ceil(gridDim/2);
    gridDimR = gridDim - gridDimL + 1;

    vl = [u10(1:gridDimL); u20(1:gridDimL)];
    vr = [u10(gridDimL:gridDim); u20(gridDimL:gridDim)];

    MLeft = MakeAdvanceMatrixLeft(C1, A, gridDimL, MakeSBPOperators, MakeBoundariesDBC);
    MRight = MakeAdvanceMatrixRight(C2, A, gridDimR, MakeSBPOperators, MakeBoundariesDBC);
    [ItLeft, ItRight] = MakeInterfaceBoundaries(A, gridDimL, gridDimR, MakeSBPOperators);

    % Do iterations
    for i=1:numIters
        ivLeft = [vl(gridDimL); vl(2*gridDimL)] - [vr(1); vr(gridDimR+1)];
        ivRight = -ivLeft;
        AdvanceLeft = @(t, vl) MLeft*vl + ItLeft*ivLeft;
        AdvanceRight = @(t, vr) MRight*vr + ItRight*ivRight;
        vl = RungeKutta4Step(AdvanceLeft, i*deltaT, vl, deltaT);
        vr = RungeKutta4Step(AdvanceRight, i*deltaT, vr, deltaT);

%         uv = [vl(1:gridDimL); vr(2:gridDimR); vl(gridDimL+1:end); vr(gridDimR+2:end)];
%         figure(9);
%         plot(x, uv(1:gridDim), 'b', x, uv(gridDim+1:end), '--r');
%         ylim([-1, 1]);
%         title(['T=', num2str(i*deltaT)]);
%         drawnow;
    end
end