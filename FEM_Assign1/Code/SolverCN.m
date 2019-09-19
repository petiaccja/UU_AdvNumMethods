function [ U, M ] = SolverCN(p, t, xiInitial, timeStep, numIters)
    C = ConvectionMatrixGFEM(p,t);
    M = MassMatrixGFEM(p,t);
    A = C/2 + M/timeStep;
    B = M/timeStep - C/2;
    
    numNodes = size(p,2);
    xi = xiInitial;

    U = zeros(numNodes, numIters+1);
    U(:,1) = xi;
    

    for i=1:numIters    
        bb = B*xi;
        xi_new = A\bb;
        xi = xi_new;
        U(:,i+1)=xi;
    end
end