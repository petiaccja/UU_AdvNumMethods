function SAT = MakeBoundariesDBC(A, H, e_1, e_m)
    tau_l = [0; 1];
    tau_r = [0; -1];
    SAT_L = kron(tau_l, inv(H)*e_1)*kron([1 0], e_1');
    SAT_R = kron(tau_r, inv(H)*e_m)*kron([1 0], e_m');
    SAT = SAT_L + SAT_R;
end