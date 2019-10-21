function [u1, u2] = AnalyticSolutionInterface(x, endT, rr, refractiveIndex1, refractiveIndex2)
    c1 = 1/refractiveIndex1;
    c2 = 1/refractiveIndex2;
    
    xWaveOriginal = -1 + (c1*endT - 0.5);
    xWaveReflected = 0 - (c1*endT - 0.5);
    xWaveTransmitted = 0 + (endT - 0.5/c1)*c2;

    Tanalytic = abs(2*refractiveIndex1/(refractiveIndex1+refractiveIndex2));
    Ranalytic = abs((refractiveIndex1-refractiveIndex2)/(refractiveIndex1+refractiveIndex2));
    
    [u1o, u2o] = InitialData(x-xWaveOriginal, rr);
    [u1t, u2t] = InitialData(x-xWaveTransmitted, rr*c2/c1);
    [u1r, u2r] = InitialData(x-xWaveReflected, rr);
    u1 = -0.5*u2o*c1 + 0.5*Ranalytic*u2r*c1 + -0.5*Tanalytic*u2t*c1;
    u2 = 0.5*u2o + 0.5*Ranalytic*u2r + 0.5*Tanalytic*u2t*c1;
    
end