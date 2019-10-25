function [ v_new ] = SAT_iterate_AD( t, v, D1, D2, DI, Hinv, eps )

    tau_l = -1;
    tau_r = 1;

    derivative = D1*v;
    d_11 = derivative(1);%[derivative(1);zeros(length(v)-1,1)];
    d_1m = derivative(end);%[zeros(length(v)-1,1);derivative(end)];
    
    e_1 = [1; zeros(length(v)-1,1)];
    e_m = [zeros(length(v)-1,1);1];

    SAT_l = tau_l*Hinv*e_1*(1/3*(v(1)+abs(v(1)))*v(1) - eps*d_11);
    SAT_r = tau_r*Hinv*e_m*(1/3*(v(end)-abs(v(end)))*v(end)-eps*d_1m);

%     v1 = v(1);
%     ve = v(end);
%     c1 = D1(1,1);
%     ce = D1(end,end);
%     SAT_l = tau_l*Hinv*e_1*(2/3*(v1 + abs(v1)) - eps*c1*v1);
%     SAT_r = tau_r*Hinv*e_m*(2/3*(ve - abs(ve)) - eps*ce*ve);
    
    %SAT_l = tau_l*Hinv * (-eps*d_11+1/3*(e_1*e_1')*v) * (e_1'*v);
    %SAT_r = tau_r*Hinv * (-eps*d_1m+1/3*(e_m*e_m')*v) * (e_m'*v);

    v_new = -1/3*D1*(v.*v) - 1/3*v.*derivative + D2*v + SAT_l + SAT_r + DI*v;

end
