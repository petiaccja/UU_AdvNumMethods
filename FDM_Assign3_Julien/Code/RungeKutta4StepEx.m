function [ ynext ] = RungeKutta4Step( f, tn, yn, ynprev, h )
    k_1 = h*f(tn, yn, ynprev);
    k_2 = h*f(tn+h/2, yn+k_1/2, ynprev);
    k_3 = h*f(tn+h/2, yn+k_2/2, ynprev);
    k_4 = h*f(tn+h, yn+k_3, ynprev);
    
    ynext = yn + 1/6*(k_1+2*k_2+2*k_3+k_4);
end
    