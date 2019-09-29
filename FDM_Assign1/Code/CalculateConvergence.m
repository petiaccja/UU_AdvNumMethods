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
        [u1, u2] = AnalyticSolution(x, endT, x_r-x_l, 0.1);
        norms(i) = sqrt(deltaX)*norm(v - [u1; u2]);
    end
    
    q = log10(norms(1)/norms(2)) / log10(gridDims(1)/gridDims(2));
end