function q = CalculateConvergenceInterface(C1, C2, A, MakeSBPOperators, MakeBoundaries)
    gridDims = [101, 201];
    x_l = -1;
    x_r = 1;
    norms = [0, 0];
    
    refractiveIndex1 = sqrt(C1(1,1));
    refractiveIndex2 = sqrt(C2(1,1));
    
    for i = 1:2
        gridDim = gridDims(i);
        deltaX = 1.0/gridDim;
        deltaT = 0.1*deltaX;
        numIters = ceil(0.42/deltaT);
        endT = deltaT*numIters;
        x = linspace(x_l, x_r, gridDim)';
        
        [vl, vr] = RunSimulationInterface(C1, C2, A, gridDim, deltaT, numIters, x_l, x_r, MakeSBPOperators, MakeBoundaries);
        gridDimL = length(vl)/2;
        gridDimR = length(vr)/2;
        v = [vl(1:gridDimL); vr(2:gridDimR); vl(gridDimL+1:end); vr(gridDimR+2:end)];
        [u1, u2] = AnalyticSolutionInterface(x, endT, 0.1, refractiveIndex1, refractiveIndex2);
        norms(i) = sqrt(deltaX)*norm(v - [u1; u2]);
        
        figure;
        plot(x, u1-v(1:gridDim), 'b', x, u2-v(gridDim+1:end), '--r');
        title(['Error: ', func2str(MakeSBPOperators), ' gridDim=', num2str(gridDim)]);
    end
    
    q = log10(norms(1)/norms(2)) / log10(gridDims(1)/gridDims(2));
end