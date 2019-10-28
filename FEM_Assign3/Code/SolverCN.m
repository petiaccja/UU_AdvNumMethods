function U = SolverCN(massMatrix, p, t, xiInitial, timeStep, numIters)
    xi = xiInitial;
    
    numNodes = length(xiInitial);
    U = zeros(numNodes, numIters+1);
    U(:,1) = xi;
    
    M = massMatrix;
    
    for i=1:numIters
        C = ConvectionMatrixGFEM(p,t,xi);
        deltaXi = (M/timeStep - C/2)\b;
        
        xi_new = xi + deltaXi;
        xi = xi_new;
        U(:,i+1)=xi;
    end
end