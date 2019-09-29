function SAT = MakeBoundariesCBC(A, H, e_1, e_m)
    [Aplus, Aminus] = SplitMatrix(A);
    SAT_R = -kron(Aplus,inv(H)*e_m)*kron(eye(2),e_m');
    SAT_L = kron(Aminus,inv(H)*e_1)*kron(eye(2),e_1');
    SAT = SAT_L + SAT_R;
end