function U = SolverCN(massMatrix, convectionMatrix, xiInitial, timeStep, numIters)
    A = convectionMatrix/2 + massMatrix/timeStep;
    B = massMatrix/timeStep - convectionMatrix/2;
    
    xi = xiInitial;
    
    numNodes = length(xiInitial);
    U = zeros(numNodes, numIters+1);
    U(:,1) = xi;    
    
    for i=1:numIters
        bb = B*xi;
        xi_new = A\bb;
        xi = xi_new;
        U(:,i+1)=xi;
    end
end