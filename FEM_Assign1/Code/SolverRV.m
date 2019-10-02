function U = SolverRV(p, t, xiInitial, timeStep, numIters)
    C = convectionMatrix;
    M = massMatrix;
    
    xi = xiInitial;
    
    numNodes = length(xiInitial);
    U = zeros(numNodes, numIters+1);
    U(:,1) = xi;
    
    r = zeros(size(t,2), 1);
    
    for i=1:numIters
        RV = ConvectionMatrixRV(p, t, xi, r);
                
        bb = B*xi;
        xi_new = A\bb;
        r = 1/deltaT*(xi_new - xi) + 
        xi = xi_new;
        U(:,i+1)=xi;
    end
end

function r = Residual(p, t, u)
    np=size(p,2); nt=size(t,2); 
    
    for i=1:nt 
        nodeIndices=t(1:3,i);
        x=p(1,nodeIndices);
        y=p(2,nodeIndices);
        z=u(nodeIndices);
        v1 = [x(1), y(1), z(1)] - [x(2), y(2), z(2)];
        
    end

end