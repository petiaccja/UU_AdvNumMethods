function convectionMatrix = ConvectionMatrixGFEM(p,t)
    numNodes=size(p,2);
    numTriangles=size(t,2);
    convectionMatrix=sparse(numNodes,numNodes);
    
    for i=1:numTriangles
        nodeIndices = t(1:3,i);
        x=p(1,nodeIndices);
        y=p(2,nodeIndices);

        [area,b,c] = HatGradients(x,y);
    
        bxmid = mean(-2*pi*y);
        bymid = mean(2*pi*x);
    
        subMatrix=ones(3,1)*(bxmid*b+bymid*c)'*area/3;
        convectionMatrix(nodeIndices,nodeIndices) = convectionMatrix(nodeIndices,nodeIndices) + subMatrix;
    end
end


function [area,b,c] = HatGradients(x,y)
    area = polyarea(x,y);
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end