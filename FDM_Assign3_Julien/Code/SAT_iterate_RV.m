function [ v_new ] = SAT_iterate_RV( t, v, vprev, D1, D2, MakeD2, Hinv, deltaT, deltaX, eps )

    tau_l = -0.5;
    tau_r = 1;

    derivative = D1*v;
    d_11 = derivative(1);
    d_1m = derivative(end);
    
    e_1 = [1; zeros(length(v)-1,1)];
    e_m = [zeros(length(v)-1,1);1];

    SAT_l = tau_l*Hinv*e_1*(1/3*(v(1)+abs(v(1)))*v(1) - eps*d_11);
    SAT_r = tau_r*Hinv*e_m*(1/3*(v(end)-abs(v(end)))*v(end)-eps*d_1m);
    
    R = Residual(v, vprev, deltaT, D1, D2);
    epsRv = EpsRv(deltaX, R, v);
    %figure(4);
    %plot(epsRv);
    %drawnow;
    
    RV = MakeD2(epsRv)*v;

    v_new = -1/3*D1*(v.*v) - 1/3*v.*derivative + D2*v + SAT_l + SAT_r + RV;
end

function R = Residual(v, vprev, deltaT, D1, D2)
    derivative = D1*v;
    R = 1/deltaT*(v-vprev) +1/3*D1*(v.*v) + 1/3*v.*derivative-D2*v;%+ 1/2*D2*(v.*v)
end

function epsRv = EpsRv(deltaX, R, v)
    Cvel = 10;
    Crv = 10;
    epsRv = min(Cvel*deltaX, Crv*deltaX^2*abs(R)/max(abs(v-mean(v))));
end
