function b = LoadAssembler2D(p,t,uprev) 
    np = size(p,2); nt = size(t,2); 
    b = zeros(np,1); 
    for K = 1:nt 
        nodeIndices = t(1:3,K); 
        x = p(1,nodeIndices);
        y = p(2,nodeIndices);
        area = polyarea(x,y); 

        fPrimePrev1 = cos(uprev(nodeIndices));
        fPrimePrev2 = -sin(uprev(nodeIndices));

        bK = (fPrimePrev1 + fPrimePrev2)/3*area; % element load vector 
        b(nodeIndices) = b(nodeIndices) + bK; % add element loads to b #
    end
end