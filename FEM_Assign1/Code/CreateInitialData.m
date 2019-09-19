function xi = CreateInitialData(p, uFunc)
    numNodes = size(p,2);
    xi = zeros(numNodes,1);
        
    for i=1:numNodes
        xi(i) = uFunc(p(1,i), p(2,i));
    end
end