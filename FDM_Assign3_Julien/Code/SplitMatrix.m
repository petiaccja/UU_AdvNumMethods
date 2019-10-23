function [Aplus, Aminus] = SplitMatrix(A) 
    [V,D] = eig(A);
    Aplus = inv(V)*(D+abs(D))/2*V;
    Aminus = inv(V)*(D-abs(D))/2*V;
end