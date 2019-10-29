function convectionMatrix = ConvectionMatrixGFEM(p,t, uprev)
    numNodes = size(p,2);
    numTriangles = size(t,2);
    convectionMatrix = sparse(numNodes, numNodes);
    
    for i=1:numTriangles
        nodeIndices = t(1:3,i);
        x=p(1,nodeIndices)';
        y=p(2,nodeIndices)';
        
        area = polyarea(x, y);
        
        fprimeUprev1 = mean(-sin(uprev(nodeIndices)));
        fprimeUprev2 = mean(-cos(uprev(nodeIndices)));
    
        subMatrix=ones(3,1)*(fprimeUprev1 + fprimeUprev2)'*area/3;
        convectionMatrix(nodeIndices,nodeIndices) = convectionMatrix(nodeIndices,nodeIndices) + subMatrix;
    end
end
