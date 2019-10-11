function [u, v] = RunSimulationInterface(C, A, gridDim, deltaT, endT, x_l, x_r, MakeSBPOperators, MakeBoundariesDBC)
    rr = 0.1;
    x = linspace(x_l, x_r, gridDim)';
    [u10, u20] = InitialData(x, rr);
    u = [u10; u20];
    v = zeros(size(u));

    MLeft = MakeAdvanceMatrixLeft(C, A, gridDim, MakeSBPOperators, MakeBoundariesDBC);
    MRight = MakeAdvanceMatrixRight(C, A, gridDim, MakeSBPOperators, MakeBoundariesDBC);
    [ItLeft, ItRight] = MakeInterfaceBoundaries(A, MakeSBPOperators);

    % Do iterations
    for i=1:ceil(endT/deltaT)
        
        AdvanceLeft = @(t, u) (M + ItLeft*1)*u;
        AdvanceRight = @(t, v) (M + ItRight*1)*v;
        u = RungeKutta4Step(AdvanceLeft, i*deltaT, u, deltaT);
        v = RungeKutta4Step(AdvanceRight, i*deltaT, v, deltaT);
    end
end