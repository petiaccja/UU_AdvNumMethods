function convectionMatrix = ConvectionMatrixRV(p, t, u, r)
    np=size(p,2); nt=size(t,2); 
    convectionMatrix=sparse(np,np);
    normFac = NormalizationFactor(p, t, u);
    for i=1:nt 
        nodeIndices=t(1:3,i); 
        x=p(1,nodeIndices); y=p(2,nodeIndices); 
        [area,b,c]=HatGradients(x,y);
        
        hk = ElementMeshSize(x, y);
        betak = ElementBetaK(x, y);
        rk = max(r(nodeIndices));
        ek = MeshArtificialViscosity(hk, betak, rk, normFac);
                
        subMatrix=(b*b')*(c+c')*area*ek;
        
        convectionMatrix(nodeIndices,nodeIndices)=convectionMatrix(nodeIndices,nodeIndices)+subMatrix;
    end
end

function hk = ElementMeshSize(x, y)
    sideXs = x - x([2,3,1]);
    sideYs = y - y([2,3,1]);
    sides = sqrt(sideXs.^2 + sideYs.^2);
    hk = max(sides);
end

function betak = ElementBetaK(x, y)
    % f'(u) = 2*pi*[-x2, x1]
    length = sqrt(x.^2 + y.^2);
    betak = 2*pi*max(length);
end

function ek = MeshArtificialViscosity(hk, betak, rk, normFac)
    Cvel = 0.25;
    Crv = 1;
    ek = min(Cvel*hk*betak, Crv*hk^2*rk/normFac);
end

function normFac = NormalizationFactor(p, t, u)
    % for calculating the mean we assume mesh elements are roughly the same
    % size
    uMean = mean(u);
    normFac = max(abs(u - uMean));
end
