function [IT_L, IT_R] = MakeInterfaceBoundaries(A, gridDimL, gridDimR, MakeSBPOperators)
    sigma_l = 1;
    sigma_r = -1;

    [Aplus, Aminus] = SplitMatrix(A);
    
    % left side
    [H, Q, D1, D2, M] = MakeSBPOperators(gridDimL);

    e_1=zeros(gridDimL,1);e_1(1)=1;
    e_m=zeros(gridDimL,1);e_m(gridDimL)=1;

    IT_L = sigma_r*kron(Aplus,inv(H)*e_m);

    % right side
    [H, Q, D1, D2, M] = MakeSBPOperators(gridDimR);

    e_1=zeros(gridDimR,1);e_1(1)=1;
    e_m=zeros(gridDimR,1);e_m(gridDimR)=1;

    IT_R = sigma_l*kron(Aminus,inv(H)*e_1);
end

