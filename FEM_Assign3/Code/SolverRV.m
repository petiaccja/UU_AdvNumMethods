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
        r = Residual(p, t, xi_new, xi, timeStep);
        xi = xi_new;
        U(:,i+1)=xi;
    end
end

function r = Residual(p, t, xi, xi_old, timeStep)
    np=size(p,2); nt=size(t,2);     
    r = ones(np, 1)*-5000;
    
    for i=1:nt 
        nodeIndices=t(1:3,i);
        
        % Calculate grad(U) for triangle.
        x=p(1,nodeIndices);
        y=p(2,nodeIndices);
        z=xi(nodeIndices);
        v1 = [x(1), y(1), z(1)] - [x(2), y(2), z(2)];
        v2 = [x(1), y(1), z(1)] - [x(3), y(3), z(3)];
        num1 = v1(3)*v2(2) - v1(2)*v2(3);
        num2 = v1(1)*v2(3) - v2(1)*v1(3);
        den = v2(2)*v1(1) - v1(2)*v2(1);
        dzdx = num1/den;
        dzdy = num2/den;
        
        % Calculate residual for the 3 points of the triangle using grad(U)
        % for the triangle.        
        z_old=xi_old(nodeIndices);
        rTriangle = 1/timeStep*(z-z_old) + 2*pi*[-y', x']*[dzdx; dzdy];
        r(nodeIndices) = max(r(nodeIndices), max(abs(rTriangle)));
    end
end