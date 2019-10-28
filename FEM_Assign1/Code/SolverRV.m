function [U, M] = SolverRV(p, t, xiInitial, timeStep, numIters)
    C = ConvectionMatrixGFEM(p, t);
    M = MassMatrixGFEM(p, t);
    
    xi = xiInitial;
    
    numNodes = length(xiInitial);
    U = zeros(numNodes, numIters+1);
    U(:,1) = xi;
    
    r = zeros(size(t,2), 1);
    
    for i=1:numIters
        RV = ConvectionMatrixRV(p, t, xi, r);
        A = 1/2*C + 1/timeStep*M + 1/2*RV;
        B = -1/2*C + 1/timeStep*M - 1/2*RV;
        bb = B*xi;
        xi_new = A\bb;
        r = Residual(M, C, xi_new, xi, timeStep);
        xi = xi_new;
        U(:,i+1)=xi;
    end
end

function r = Residual(M, C, xi, xi_old, timeStep)
    lhs = 1/timeStep*(M*(xi - xi_old)) + C*xi;
    r = M\lhs;
end