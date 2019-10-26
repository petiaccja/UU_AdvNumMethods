function convectionMatrix = ConvectionMatrixGFEM(p,t, uprev)
    numNodes = size(p,2);
    numTriangles = size(t,2);
    convectionMatrix = sparse(numNodes, numNodes);
    
    for i=1:numTriangles
        nodeIndices = t(1:3,i);
        x=p(1,nodeIndices);
        y=p(2,nodeIndices);

        [area,b,c] = HatGradients(x,y);
        [sK, cK] = fK(nodeIndices, uprev);
            
        subMatrix=ones(3,1)*(sK*b + cK*c)'*area/3;
        convectionMatrix(nodeIndices,nodeIndices) = convectionMatrix(nodeIndices,nodeIndices) + subMatrix;
    end
end

function [s,c] = f(u)
    s = sin(u);
    c = cos(u);
end

function [sK, cK] = fK(nodeIndices, u)
    [s, c] = f(u(nodeIndices));
    sK = mean(s);
    cK = mean(c);
end
