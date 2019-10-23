function M = MakeAdvanceMatrix(C, A, gridDim, MakeSBPOperators, MakeBoundaries)
    [H, Q, D1, D2, M] = MakeSBPOperators(gridDim);

    e_1=zeros(gridDim,1);e_1(1)=1;
    e_m=zeros(gridDim,1);e_m(gridDim)=1;
    
    Cmod = kron(C, eye(gridDim));
    FD = kron(A, D1);
    SAT = MakeBoundaries(A, H, e_1, e_m);
    M = inv(Cmod)*(FD+SAT);
end