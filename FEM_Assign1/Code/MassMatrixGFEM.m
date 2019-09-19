function massMatrix = MassMatrixGFEM(p,t)
    numNodes = size(p,2);
    numTriangles = size(t,2);
    massMatrix = sparse(numNodes,numNodes);

    for k = 1:numTriangles
        nodeIndices = t(1:3, k);
        x = p(1,nodeIndices); % Node x-coordinates
        y = p(2,nodeIndices); % Node y
        triangleArea = polyarea(x,y);
        subMatrix = [2,1,1;1,2,1;1,1,2]/12*triangleArea;
        massMatrix(nodeIndices,nodeIndices) = massMatrix(nodeIndices,nodeIndices) + subMatrix;
    end
end