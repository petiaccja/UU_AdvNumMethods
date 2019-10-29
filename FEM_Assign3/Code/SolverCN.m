function U = SolverCN(massMatrix, p, t, xiInitial, timeStep, numIters)
    xi = xiInitial;
    
    numNodes = length(xiInitial);
    U = zeros(numNodes, numIters+1);
    U(:,1) = xi;
    
    M = massMatrix;
    
    for i=1:numIters
        C = ConvectionMatrixGFEM(p,t,xi);
        b = LoadAssembler2D(p,t,xi);
        
        Zeta = zeros(size(xi));
        
        for j=1:20
            Zeta = -timeStep*M\(b+1/2*C*Zeta);
        end %for, j
        
        deltaXi = Zeta;
        
        deltaXi_questionable = (M/timeStep - C/2)\b;
        
        %figure(6);
        %pdesurf(p, t, deltaXi);
        %drawnow;
        
        xi_new = xi + deltaXi;
        xi = xi_new;
        U(:,i+1)=xi;
    end
end