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