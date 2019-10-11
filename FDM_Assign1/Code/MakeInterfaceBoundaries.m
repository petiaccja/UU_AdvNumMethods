function [IT_L, IT_R] = MakeInterfaceBoundaries(A, MakeSBPOperators)
    [H, Q, D1, D2, M] = MakeSBPOperators(gridDim);

    e_1=zeros(gridDim,1);e_1(1)=1;
    e_m=zeros(gridDim,1);e_m(gridDim)=1;

    sigma_l = 1;
    sigma_r = -1;

    [Aplus, Aminus] = SplitMatrix(A);

    IT_L = sigma_r*kron(Aplus,inv(H)*e_m);
    IT_R = sigma_l*kron(Aminus,inv(H)*e_1);
end

