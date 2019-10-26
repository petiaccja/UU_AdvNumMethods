function U = SolverCN(massMatrix, p, t, xiInitial, timeStep, numIters)
    xi = xiInitial;
    
    numNodes = length(xiInitial);
    U = zeros(numNodes, numIters+1);
    U(:,1) = xi;    
    
    for i=1:numIters
        convectionMatrix = ConvectionMatrixGFEM(p,t,xi);
        
        A = convectionMatrix/2 + massMatrix/timeStep;
        B = massMatrix/timeStep - convectionMatrix/2;
    
        bb = B*xi;
        xi_new = A\bb;
        xi = xi_new;
        U(:,i+1)=xi;
    end
end